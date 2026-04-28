# scGPT Model Setup

## Pretrained Models

8 pretrained models available, trained on different cell populations:

| Model | Training Data | Size | Source |
|-------|--------------|------|--------|
| **whole-human** (recommended) | 33M normal human cells | ~208MB | [Google Drive](https://drive.google.com/drive/folders/1oWh_-ZRdhtoGQ2Fw24HP41FgLoomVo-y) |
| continual pretrained | Zero-shot cell embedding | ~208MB | [Google Drive](https://drive.google.com/drive/folders/1_GROJTzXiAV8HB4imruOTk6PEGuNOcgB) |
| brain | 13.2M brain cells | ~208MB | [Google Drive](https://drive.google.com/drive/folders/1vf1ijfQSk7rGdDGpBntR5bi5g6gNt-Gx) |
| blood | 10.3M blood/bone marrow | ~208MB | [Google Drive](https://drive.google.com/drive/folders/1kkug5C7NjvXIwQGGaGoqXTk_Lb_pDrBU) |
| heart | 1.8M heart cells | ~208MB | [Google Drive](https://drive.google.com/drive/folders/1GcgXrd7apn6y4Ze_iSCncskX3UsWPY2r) |
| lung | 2.1M lung cells | ~208MB | [Google Drive](https://drive.google.com/drive/folders/16A1DJ30PT6bodt4bWLa4hpS7gbWZQFBG) |
| kidney | 814K kidney cells | ~208MB | [Google Drive](https://drive.google.com/drive/folders/1S-1AR65DF120kNFpEbWCvRHPhpkGK3kK) |
| pan-cancer | 5.7M cancer cells | ~208MB | [Google Drive](https://drive.google.com/drive/folders/13QzLHilYUd0v3HTwa_9n4G4yEF-hdkqa) |

## Download (Automated)

Use the download script:

```bash
eval "$(conda shell.bash hook 2>/dev/null)" && conda activate scgpt
python /home/sagemaker-user/.claude/skills/scgpt/scripts/download_model.py \
  --output-dir /path/to/destination
```

This downloads the whole-human model from a HuggingFace community mirror (`MohamedMabrouk/scGPT`).

## Download (Manual from Google Drive)

For tissue-specific models, use gdown:

```bash
pip install gdown
# whole-human model
gdown --folder "https://drive.google.com/drive/folders/1oWh_-ZRdhtoGQ2Fw24HP41FgLoomVo-y" -O /path/to/scgpt_model/
```

## Current Local Path

```
/home/sagemaker-user/claude_code/scimilarity/scgpt_model/
├── best_model.pt    # Pretrained weights (~208MB)
├── vocab.json       # Gene vocabulary (48,292 HGNC symbols)
└── args.json        # Model configuration
```

## S3 Storage & Cache Helper (Recommended)

Model is persisted on S3. Use the cache helper to auto-sync to `/tmp/scgpt/`:

```python
import sys; sys.path.insert(0, "/home/sagemaker-user/.claude/skills/scgpt/scripts")
from s3_cache import get_model_path

model_dir = get_model_path()  # syncs ~210MB on first call, instant after
# model_dir = "/tmp/scgpt"
```

Or from CLI:
```bash
python /home/sagemaker-user/.claude/skills/scgpt/scripts/s3_cache.py
python /home/sagemaker-user/.claude/skills/scgpt/scripts/s3_cache.py --clear  # remove cache
```

All pipeline scripts support `--s3`:
```bash
python embed_cells.py --input data.h5ad --s3 --output embedded.h5ad
python finetune.py --input data.h5ad --s3 --task annotation --output-dir ./finetuned
```

## Model Files

| File | Purpose | Size |
|------|---------|------|
| `best_model.pt` | PyTorch checkpoint with model state_dict | ~208MB |
| `vocab.json` | Gene name to integer ID mapping (48,292 genes) | ~1.3MB |
| `args.json` | Model hyperparameters (d_model, nhead, nlayers, etc.) | <1KB |

## Loading the Model

```python
import json
import torch
from scgpt.model.model import TransformerModel
from scgpt.tokenizer import GeneVocab
from scgpt.utils import load_pretrained

model_dir = "/home/sagemaker-user/claude_code/scimilarity/scgpt_model"

# Load vocab
vocab = GeneVocab.from_file(f"{model_dir}/vocab.json")

# Load args
with open(f"{model_dir}/args.json") as f:
    model_args = json.load(f)

# Build model
model = TransformerModel(
    ntoken=len(vocab),
    d_model=model_args.get("embsize", 512),
    nhead=model_args.get("nheads", 8),
    d_hid=model_args.get("d_hid", 512),
    nlayers=model_args.get("nlayers", 12),
    vocab=vocab,
    pad_token="<pad>",
    # Add task-specific flags as needed
)

# Load pretrained weights
checkpoint = torch.load(f"{model_dir}/best_model.pt", map_location="cpu")
load_pretrained(model, checkpoint)
model.eval()
```

## Verification

```python
eval "$(conda shell.bash hook 2>/dev/null)" && conda activate scgpt
python -c "
import torch, json
from scgpt.tokenizer import GeneVocab
vocab = GeneVocab.from_file('scgpt_model/vocab.json')
print(f'Vocab: {len(vocab)} genes')
ckpt = torch.load('scgpt_model/best_model.pt', map_location='cpu')
print(f'Checkpoint keys: {len(ckpt)} layers')
with open('scgpt_model/args.json') as f:
    args = json.load(f)
print(f'Model args: {args}')
"
```

## Default Architecture (whole-human / continual pretrained)

- **d_model**: 512 (embedding dimension)
- **nhead**: 8 (attention heads)
- **d_hid**: 512 (FFN hidden dim)
- **nlayers**: 12 (transformer layers)
- **n_cls**: 177 (pretrained annotation classes)
- **Parameters**: ~51.9M total
- **Vocabulary**: 60,697 tokens (genes + special: `<pad>`, `<cls>`, `<eoc>`)
- **Max sequence length**: 1200 genes per cell
- **pad_value**: -2 (NOT 0 — this model uses -2 for padding)
- **Normalization**: CPM 1e4 + log1p (same as scimilarity)
- **Input style**: Binned (51 bins) with continuous value encoder

## Important: Vocab Loading

The `GeneVocab.from_file()` may fail with this checkpoint's vocab.json due to a torchtext compatibility issue. Use a simple dict wrapper instead:

```python
import json

class SimpleVocab:
    def __init__(self, token2idx):
        self.token2idx = token2idx
    def __getitem__(self, token):
        return self.token2idx[token]
    def __len__(self):
        return len(self.token2idx)
    def __contains__(self, token):
        return token in self.token2idx

with open("scgpt_model/vocab.json") as f:
    vocab = SimpleVocab(json.load(f))
```
