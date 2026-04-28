#!/usr/bin/env python3
"""CLI wrapper for Geneformer MTLClassifier (multi-task learning).

Fine-tunes Geneformer on multiple classification tasks simultaneously
with Optuna hyperparameter optimization.

Usage:
    conda activate geneformer

    # Train with Optuna hyperparameter search
    python mtl_classify.py train --s3 \
        --task-columns cell_type,disease \
        --train-path train.dataset --val-path val.dataset \
        --model-save-path mtl_model/ --results-dir mtl_results/ \
        --study-name my_mtl --n-trials 15 --epochs 10

    # Train with attention pooling and task weights
    python mtl_classify.py train --s3 --s3-tier V2-316M \
        --task-columns cell_type,disease,organ \
        --train-path train.dataset --val-path val.dataset \
        --model-save-path mtl_model/ --results-dir mtl_results/ \
        --study-name my_mtl \
        --use-attention-pooling --use-task-weights

    # Evaluate a trained model on test data
    python mtl_classify.py evaluate \
        --task-columns cell_type,disease \
        --test-path test.dataset \
        --model-save-path mtl_model/ --results-dir mtl_results/

    # Train with direct model path
    python mtl_classify.py train \
        --pretrained-path /path/to/pretrained/model \
        --task-columns cell_type,disease \
        --train-path train.dataset --val-path val.dataset \
        --model-save-path mtl_model/ --results-dir mtl_results/ \
        --study-name my_mtl
"""

import argparse
import os
import sys


def parse_args():
    parser = argparse.ArgumentParser(
        description="Geneformer multi-task classification with Optuna",
        formatter_class=argparse.RawDescriptionHelpFormatter,
    )

    # Subcommand
    subparsers = parser.add_subparsers(dest="command", required=True)

    # Shared arguments for both subcommands
    shared = argparse.ArgumentParser(add_help=False)
    shared.add_argument(
        "--task-columns",
        required=True,
        type=str,
        help="Comma-separated task column names from the dataset",
    )
    shared.add_argument("--model-save-path", required=True, help="Directory to save/load the fine-tuned model")
    shared.add_argument("--results-dir", required=True, help="Directory for results and trial logs")
    shared.add_argument("--seed", type=int, default=42, help="Random seed (default: 42)")

    # Train subcommand
    train_parser = subparsers.add_parser("train", parents=[shared], help="Train with Optuna hyperparameter search")

    # Model source for training (mutually exclusive)
    model_group = train_parser.add_mutually_exclusive_group(required=True)
    model_group.add_argument("--pretrained-path", type=str, help="Direct path to pretrained model")
    model_group.add_argument("--s3", action="store_true", help="Auto-cache pretrained model from S3")
    train_parser.add_argument(
        "--s3-tier",
        type=str,
        default="V2-104M",
        help="S3 model tier (default: V2-104M)",
    )

    # Training data
    train_parser.add_argument("--train-path", required=True, help="Path to training .dataset")
    train_parser.add_argument("--val-path", required=True, help="Path to validation .dataset")

    # Training hyperparameters
    train_parser.add_argument("--study-name", type=str, default=None, help="Optuna study name")
    train_parser.add_argument("--batch-size", type=int, default=4, help="Training batch size (default: 4)")
    train_parser.add_argument("--n-trials", type=int, default=15, help="Number of Optuna trials (default: 15)")
    train_parser.add_argument("--epochs", type=int, default=None, help="Number of training epochs")
    train_parser.add_argument(
        "--use-attention-pooling",
        action="store_true",
        help="Use attention-based pooling (recommended)",
    )
    train_parser.add_argument(
        "--use-task-weights",
        action="store_true",
        help="Use learnable task weights",
    )
    train_parser.add_argument(
        "--distributed-training",
        action="store_true",
        help="Enable distributed data parallel training",
    )
    train_parser.add_argument(
        "--gradient-accumulation-steps",
        type=int,
        default=None,
        help="Gradient accumulation steps",
    )

    # Evaluate subcommand
    eval_parser = subparsers.add_parser("evaluate", parents=[shared], help="Evaluate a trained model on test data")
    eval_parser.add_argument("--test-path", required=True, help="Path to test .dataset")

    return parser.parse_args()


def main():
    args = parse_args()

    # Parse task columns
    task_columns = [col.strip() for col in args.task_columns.split(",")]

    if args.command == "train":
        # Resolve pretrained model path
        if args.s3:
            sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
            from s3_cache import get_model_path

            pretrained_path = get_model_path(tier=args.s3_tier)
        else:
            pretrained_path = args.pretrained_path

        # Build trials result path
        trials_result_path = os.path.join(args.results_dir, "trials.txt")
        os.makedirs(args.results_dir, exist_ok=True)

        from geneformer import MTLClassifier

        mc = MTLClassifier(
            task_columns=task_columns,
            pretrained_path=pretrained_path,
            train_path=args.train_path,
            val_path=args.val_path,
            model_save_path=args.model_save_path,
            results_dir=args.results_dir,
            trials_result_path=trials_result_path,
            study_name=args.study_name,
            batch_size=args.batch_size,
            n_trials=args.n_trials,
            epochs=args.epochs,
            use_attention_pooling=args.use_attention_pooling or None,
            use_task_weights=args.use_task_weights or None,
            distributed_training=args.distributed_training or None,
            seed=args.seed,
            gradient_accumulation_steps=args.gradient_accumulation_steps,
        )

        print(f"Starting Optuna study with {args.n_trials} trials...")
        print(f"Tasks: {task_columns}")
        print(f"Pretrained model: {pretrained_path}")
        mc.run_optuna_study()
        print(f"Training complete. Model saved to: {args.model_save_path}")
        print(f"Results saved to: {args.results_dir}")

    elif args.command == "evaluate":
        from geneformer import MTLClassifier

        os.makedirs(args.results_dir, exist_ok=True)

        mc = MTLClassifier(
            task_columns=task_columns,
            test_path=args.test_path,
            model_save_path=args.model_save_path,
            results_dir=args.results_dir,
            seed=args.seed,
        )

        print(f"Evaluating model from: {args.model_save_path}")
        print(f"Test data: {args.test_path}")
        mc.load_and_evaluate_test_model()
        print(f"Evaluation complete. Results saved to: {args.results_dir}")


if __name__ == "__main__":
    main()
