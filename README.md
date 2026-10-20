# Nucleo — DNA Sequence Visualizer

A Python web app for introductory DNA sequence analysis.

## Features
1 Complementary DNA strand
2 Reverse complement
3 DNA-to-RNA transcription (coding-strand convention)
4 Translation using the standard genetic code (one- and three-letter amino-acid codes)
5 GC content and A/T/G/C composition chart
6 Reading-frame selection and optional stop-codon termination
7 Downloadable TXT report

## Run locally

Python 3.10+ recommended.

```bash
python -m venv .venv
# macOS / Linux
source .venv/bin/activate
# Windows PowerShell
# .venv\Scripts\Activate.ps1

pip install -r requirements.txt
streamlit run app.py
```

The app opens in your browser, usually at `http://localhost:8501`.

## Notes
1 Input accepts only A, T, C, G; whitespace is ignored.
2 RNA conversion assumes the input is the coding DNA strand, so it replaces T with U.
3 Translation uses the standard genetic code. It translates the selected reading frame and does not automatically search for a start codon.
4 This is an educational tool, not a clinical or diagnostic application.
