# Nucleo — DNA Sequence Visualizer

A dark-mode Python web app for introductory DNA sequence analysis.

## Features
- Complementary DNA strand
- Reverse complement
- DNA-to-RNA transcription (coding-strand convention)
- Translation using the standard genetic code (one- and three-letter amino-acid codes)
- GC content and A/T/G/C composition chart
- Reading-frame selection and optional stop-codon termination
- Downloadable TXT report

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
- Input accepts only A, T, C, G; whitespace is ignored.
- RNA conversion assumes the input is the coding DNA strand, so it replaces T with U.
- Translation uses the standard genetic code. It translates the selected reading frame and does not automatically search for a start codon.
- This is an educational tool, not a clinical or diagnostic application.
