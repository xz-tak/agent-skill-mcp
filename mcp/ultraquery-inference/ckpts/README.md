# ULTRA Model Checkpoints

This directory stores pre-trained ULTRA model checkpoints for inference.

## Automatic Download

**Checkpoints are automatically downloaded from GitHub** when you first run inference. No manual setup needed!

The MCP server will:
1. Check if the checkpoint (`ultraquery_primekg_50g_ft_epoch_1.pth`) exists
2. Download it from GitHub if not present
3. Verify the download was successful
4. Load the model and begin inference

Progress is logged during download so you can track the status.

## Checkpoint Details

The checkpoint (`ultraquery_primekg_50g_ft_epoch_1.pth`) is:
- **Base model**: ULTRAQuery 50g (pre-trained on 50 knowledge graphs)
- **Fine-tuned on**: PrimeKG biomedical knowledge graph
- **Size**: ~200-500 MB (large file, excluded from git)
- **Architecture**: Dual-graph NBFNet with RelNBFNet + EntityNBFNet
- **Source**: https://github.com/roger-tu/ULTRA

## Note

Checkpoint files (*.pth) are excluded from version control via `.gitignore` due to their large size. The automatic download ensures the checkpoint is always available when needed.
