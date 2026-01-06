#!/usr/bin/env python
from __future__ import annotations

import argparse
import pandas as pd


def main():
    ap = argparse.ArgumentParser(description='Inspect an Excel workbook (sheets, columns, previews).')
    ap.add_argument('xlsx', help='Path to .xlsx file')
    ap.add_argument('--rows', type=int, default=3, help='Preview rows per sheet')
    args = ap.parse_args()

    xl = pd.ExcelFile(args.xlsx)
    print('Sheets:', xl.sheet_names)
    for name in xl.sheet_names:
        df = pd.read_excel(args.xlsx, sheet_name=name)
        print(f"\n== {name} ==")
        print('shape:', df.shape)
        print('columns:', list(df.columns))
        if len(df):
            print(df.head(args.rows).to_string(index=False))


if __name__ == '__main__':
    main()
