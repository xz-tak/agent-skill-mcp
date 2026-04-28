#!/usr/bin/env python3
"""Geneformer comprehensive workflow examples.

Each function demonstrates a self-contained Geneformer workflow with
placeholder paths. Run any example from the command line:

Usage:
    conda activate geneformer

    python comprehensive_workflow.py tokenize
    python comprehensive_workflow.py embed
    python comprehensive_workflow.py classify
    python comprehensive_workflow.py perturb
    python comprehensive_workflow.py mtl
    python comprehensive_workflow.py quantize
"""

import argparse
import os
import sys

# Add scripts dir to path for s3_cache
SCRIPTS_DIR = os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "scripts")
sys.path.insert(0, SCRIPTS_DIR)

from s3_cache import get_model_path


def tokenize_example():
    """Step 1: Tokenize raw scRNA-seq h5ad files into rank value encodings.

    Input: directory of .h5ad files with ensembl_id and n_counts
    Output: HuggingFace .dataset with input_ids, length, and metadata
    """
    from geneformer import TranscriptomeTokenizer

    # -- Configure tokenizer --
    # Map h5ad obs columns to dataset column names
    custom_attrs = {
        "cell_type": "cell_type",   # cell type annotation column
        "organ_major": "organ",     # organ/tissue column
    }

    tk = TranscriptomeTokenizer(
        custom_attr_name_dict=custom_attrs,
        nproc=16,
        chunk_size=512,
        model_input_size=4096,      # V2 default, do not change
        special_token=True,         # V2 default, adds CLS/EOS tokens
        collapse_gene_ids=True,     # merge duplicate Ensembl IDs
    )

    # -- Tokenize all h5ad files in directory --
    # Replace with your actual paths
    data_directory = "/path/to/h5ad/files"       # directory containing .h5ad files
    output_directory = "/path/to/output"          # where tokenized .dataset will be saved
    output_prefix = "my_tokenized"                # prefix for output files

    tk.tokenize_data(
        data_directory,
        output_directory,
        output_prefix,
        file_format="h5ad",
    )
    print(f"Tokenized data saved to: {output_directory}/{output_prefix}.dataset")


def embed_example():
    """Step 2: Extract cell embeddings and generate UMAP visualization.

    Input: tokenized .dataset from Step 1
    Output: embedding dataframe (pandas) + UMAP plot
    """
    from geneformer import EmbExtractor

    # -- Load pretrained model from S3 cache --
    model_path = get_model_path(tier="V2-104M")

    # -- Configure embedding extractor --
    embex = EmbExtractor(
        model_type="Pretrained",
        num_classes=0,              # 0 for pretrained (no classification head)
        emb_mode="cls",             # CLS token embedding (V2 recommended)
        max_ncells=2000,            # limit for faster exploration; None for all cells
        emb_layer=-1,               # -1: second-to-last layer (general purpose)
        emb_label=["cell_type"],    # columns to include as labels in output
        labels_to_plot=["cell_type"],
        forward_batch_size=256,
        nproc=16,
    )

    # -- Extract embeddings --
    # Replace with your actual tokenized dataset path
    input_data = "/path/to/my_tokenized.dataset"
    output_directory = "/path/to/embeddings"
    output_prefix = "my_embs"

    embs = embex.extract_embs(
        model_path,
        input_data,
        output_directory,
        output_prefix,
    )
    print(f"Extracted embeddings shape: {embs.shape}")

    # -- Generate UMAP plot --
    embex.plot_embs(
        embs,
        plot_style="umap",
        output_directory=output_directory,
        output_prefix=output_prefix,
    )
    print(f"UMAP plot saved to: {output_directory}")


def classify_example():
    """Step 3: Fine-tune a cell state classifier.

    Input: tokenized .dataset with a state column (e.g., disease)
    Output: trained model, confusion matrix, metrics
    """
    from geneformer import Classifier

    # -- Load pretrained model from S3 cache --
    model_path = get_model_path(tier="V2-104M")

    # -- Configure classifier --
    cc = Classifier(
        classifier="cell",          # "cell" for cell-level, "gene" for gene-level
        cell_state_dict={
            "state_key": "disease",                   # obs column to classify
            "states": ["healthy", "diseased"],        # classes to distinguish
        },
        filter_data={"cell_type": ["Cardiomyocyte"]}, # subset to specific cells
        training_args={
            "num_train_epochs": 3,
            "learning_rate": 5e-5,
            "per_device_train_batch_size": 12,
            "warmup_steps": 500,
        },
        freeze_layers=2,            # freeze first 2 transformer layers
        num_crossval_splits=1,      # 1: train/eval split, 5: 5-fold CV
        forward_batch_size=200,
        nproc=16,
    )

    # -- Prepare data (label and split) --
    input_data = "/path/to/my_tokenized.dataset"
    output_directory = "/path/to/classifier_output"
    output_prefix = "my_classifier"

    cc.prepare_data(input_data, output_directory, output_prefix)

    # -- Train and cross-validate --
    labeled_data = os.path.join(output_directory, f"{output_prefix}_labeled.dataset")
    id_class_dict = os.path.join(output_directory, f"{output_prefix}_id_class_dict.pkl")

    all_metrics = cc.validate(
        model_path,
        labeled_data,
        id_class_dict,
        output_directory,
        output_prefix,
        predict_eval=True,
    )

    # -- Plot confusion matrix --
    cc.plot_conf_mat(
        {"Geneformer": all_metrics["conf_matrix"]},
        output_directory,
        output_prefix,
    )
    print(f"Classifier results saved to: {output_directory}")


def perturb_example():
    """Step 4: Full in silico perturbation pipeline (3 steps).

    Simulates deleting every gene and measures which deletions shift
    diseased cells toward the healthy (goal) state.

    Input: tokenized .dataset, fine-tuned model
    Output: ranked gene list (CSV) by goal state shift
    """
    from geneformer import EmbExtractor, InSilicoPerturber, InSilicoPerturberStats

    # -- Load model from S3 cache --
    model_path = get_model_path(tier="V2-104M")

    # -- Define cell states for reverse phenotype analysis --
    cell_states = {
        "state_key": "disease",
        "start_state": "dcm",       # starting (disease) state
        "goal_state": "nf",          # goal (healthy) state
        "alt_states": ["hcm"],       # alternative states to track
    }

    # Replace with your actual paths
    input_data = "/path/to/my_tokenized.dataset"
    isp_output_dir = "/path/to/isp_output"
    stats_output_dir = "/path/to/stats_output"

    # ==========================================
    # Step 4a: Compute state embeddings
    # ==========================================
    embex = EmbExtractor(
        model_type="CellClassifier",
        num_classes=3,               # must match number of states (dcm, nf, hcm)
        filter_data={"cell_type": ["Cardiomyocyte"]},
        max_ncells=1000,
        emb_layer=0,                 # 0: last layer (task-specific for fine-tuned)
        summary_stat="exact_mean",   # compute exact mean for state centroids
        forward_batch_size=256,
        nproc=16,
    )

    state_embs = embex.get_state_embs(
        cell_states,
        model_path,
        input_data,
        isp_output_dir,
        "state_embs",
    )
    print(f"Computed state embeddings for {len(state_embs)} states")

    # ==========================================
    # Step 4b: Run perturbation (delete each gene)
    # ==========================================
    isp = InSilicoPerturber(
        perturb_type="delete",       # delete each gene one at a time
        genes_to_perturb="all",      # perturb all genes (or pass a list of Ensembl IDs)
        model_type="CellClassifier",
        num_classes=3,
        emb_mode="cls",              # use CLS token embedding for shift measurement
        filter_data={"cell_type": ["Cardiomyocyte"]},
        cell_states_to_model=cell_states,
        state_embs_dict=state_embs,
        max_ncells=2000,
        emb_layer=0,
        forward_batch_size=400,
        nproc=16,
    )

    isp.perturb_data(model_path, input_data, isp_output_dir, "my_isp")
    print(f"Perturbation results saved to: {isp_output_dir}")

    # ==========================================
    # Step 4c: Compute statistics (goal state shift)
    # ==========================================
    ispstats = InSilicoPerturberStats(
        mode="goal_state_shift",     # rank genes by shift toward goal state
        genes_perturbed="all",
        cell_states_to_model=cell_states,
    )

    ispstats.get_stats(isp_output_dir, None, stats_output_dir, "my_stats")
    print(f"Gene rankings saved to: {stats_output_dir}")


def mtl_example():
    """Step 5: Multi-task classification with Optuna hyperparameter search.

    Trains on multiple classification tasks simultaneously (e.g., cell_type
    and disease) with automatic hyperparameter optimization.

    Input: tokenized .dataset with unique_cell_id + task columns
    Output: fine-tuned model, Optuna trial results, test evaluation
    """
    from geneformer import MTLClassifier

    # -- Load pretrained model from S3 cache --
    model_path = get_model_path(tier="V2-104M")

    # Replace with your actual paths
    train_path = "/path/to/train.dataset"
    val_path = "/path/to/val.dataset"
    test_path = "/path/to/test.dataset"
    model_save_path = "/path/to/mtl_model"
    results_dir = "/path/to/mtl_results"

    # -- Configure and train --
    mc = MTLClassifier(
        task_columns=["cell_type", "disease"],   # dataset columns for each task
        pretrained_path=model_path,
        train_path=train_path,
        val_path=val_path,
        test_path=test_path,
        model_save_path=model_save_path,
        results_dir=results_dir,
        trials_result_path=os.path.join(results_dir, "trials.txt"),
        study_name="mtl_study",
        batch_size=8,
        n_trials=15,                 # number of Optuna hyperparameter trials
        epochs=10,
        use_attention_pooling=True,  # attention-based pooling (recommended)
        use_task_weights=True,       # learnable task weights
        seed=42,
    )

    # -- Run Optuna study (hyperparameter search + training) --
    mc.run_optuna_study()
    print(f"Training complete. Model saved to: {model_save_path}")

    # -- Evaluate best model on test data --
    mc.load_and_evaluate_test_model()
    print(f"Test evaluation saved to: {results_dir}")


def quantize_example():
    """Step 6: Fine-tune with 8-bit quantization for memory-efficient training.

    Uses BitsAndBytesConfig for 8-bit quantization and LoRA for
    parameter-efficient fine-tuning. Reduces GPU memory by ~50%.

    Input: tokenized .dataset
    Output: quantized fine-tuned model
    """
    from geneformer import Classifier

    # -- Load pretrained model from S3 cache --
    model_path = get_model_path(tier="V2-104M")

    # -- Option A: Simple 8-bit quantization --
    cc_simple = Classifier(
        classifier="cell",
        quantize=True,               # enables 8-bit quantization
        cell_state_dict={
            "state_key": "disease",
            "states": ["healthy", "diseased"],
        },
        training_args={
            "num_train_epochs": 3,
            "learning_rate": 5e-5,
            "per_device_train_batch_size": 8,
        },
        freeze_layers=2,
        num_crossval_splits=1,
        forward_batch_size=200,
        nproc=16,
    )
    print("Simple quantized classifier created (quantize=True)")

    # -- Option B: Custom BitsAndBytes + LoRA configuration --
    from peft import LoraConfig
    from transformers import BitsAndBytesConfig

    cc_custom = Classifier(
        classifier="cell",
        quantize={
            "bnb_config": BitsAndBytesConfig(load_in_8bit=True),
            "peft_config": LoraConfig(
                r=8,
                lora_alpha=32,
                target_modules=["query", "value"],
            ),
        },
        cell_state_dict={
            "state_key": "disease",
            "states": ["healthy", "diseased"],
        },
        training_args={
            "num_train_epochs": 3,
            "learning_rate": 5e-5,
            "per_device_train_batch_size": 8,
        },
        freeze_layers=0,            # LoRA handles parameter efficiency
        num_crossval_splits=1,
        forward_batch_size=200,
        nproc=16,
    )
    print("Custom quantized classifier created (BnB + LoRA)")

    # -- Prepare data and train (same API as non-quantized) --
    input_data = "/path/to/my_tokenized.dataset"
    output_directory = "/path/to/quantized_output"
    output_prefix = "quant_classifier"

    cc_simple.prepare_data(input_data, output_directory, output_prefix)

    labeled_data = os.path.join(output_directory, f"{output_prefix}_labeled.dataset")
    id_class_dict = os.path.join(output_directory, f"{output_prefix}_id_class_dict.pkl")

    all_metrics = cc_simple.validate(
        model_path,
        labeled_data,
        id_class_dict,
        output_directory,
        output_prefix,
        predict_eval=True,
    )
    print(f"Quantized classifier results saved to: {output_directory}")


def main():
    parser = argparse.ArgumentParser(
        description="Geneformer comprehensive workflow examples",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Examples:
    python comprehensive_workflow.py tokenize    # Tokenize raw h5ad data
    python comprehensive_workflow.py embed       # Extract cell embeddings + UMAP
    python comprehensive_workflow.py classify    # Fine-tune cell state classifier
    python comprehensive_workflow.py perturb     # In silico perturbation pipeline
    python comprehensive_workflow.py mtl         # Multi-task classification
    python comprehensive_workflow.py quantize    # Quantized fine-tuning
        """,
    )
    parser.add_argument(
        "example",
        choices=["tokenize", "embed", "classify", "perturb", "mtl", "quantize"],
        help="Which workflow example to run",
    )

    args = parser.parse_args()

    dispatch = {
        "tokenize": tokenize_example,
        "embed": embed_example,
        "classify": classify_example,
        "perturb": perturb_example,
        "mtl": mtl_example,
        "quantize": quantize_example,
    }

    dispatch[args.example]()


if __name__ == "__main__":
    main()
