#!/usr/bin/env python3
"""
Build ASI1, ASI2, ASI3 supplementary PDFs (vertical A4).

Architecture:
- Source PDFs -> PNG (300 dpi)
- Each page fully composed with ReportLab
- PyPDF2 used ONLY to concatenate finished pages
- Robust to corrupted PDFs
"""

import re
from pathlib import Path
from io import BytesIO

import pandas as pd
from tqdm import tqdm

from reportlab.lib.pagesizes import A4
from reportlab.pdfgen import canvas
from reportlab.lib.units import cm

from PyPDF2 import PdfReader, PdfWriter
from pdf2image import convert_from_path
from PIL import Image, ImageFile

ImageFile.LOAD_TRUNCATED_IMAGES = True

# =============================================================================
# PATHS
# =============================================================================
BASE = Path("/Users/jagon/Documents/Projects/Collabs/Jessica Badgeley")
SI_DIR = BASE / "Manuscript/SI"
CSV_PATH = BASE / "New Points v3/Input/Box Coordinates/box_sp_all_v3.csv"

ASI_FRONT = {
    1: SI_DIR / "ASI1_front_page.pdf",
    2: SI_DIR / "ASI2_front_page.pdf",
    3: SI_DIR / "ASI3_front_page.pdf",
}

ASI1_INPUT = BASE / "New Points v3/Flowlines of Interest/Graphs with Picks/with Lines Chad"

ASI2_INPUT = {
    "fric1": BASE / "Version 3/Outputs/fric1/Graphs/Obs-Mod-Res",
    "fric2": BASE / "Version 3/Outputs/fric2/Graphs/Obs-Mod-Res",
    "fric3": BASE / "Version 3/Outputs/fric3/Graphs/Obs-Mod-Res",
}

ASI3_INPUT = {
    "fric1": BASE / "Version 3/Outputs/fric1/Heatmaps/Regular",
    "fric2": BASE / "Version 3/Outputs/fric2/Heatmaps/Regular",
    "fric3": BASE / "Version 3/Outputs/fric3/Heatmaps/Regular",
}

OUTFILES = {
    1: SI_DIR / "ASI1.pdf",
    2: SI_DIR / "ASI2.pdf",
    3: SI_DIR / "ASI3.pdf",
}

PNG_CACHE = SI_DIR / "_png_cache"
PNG_CACHE.mkdir(exist_ok=True)

PAGE_W, PAGE_H = A4
PNG_DPI = 300

LABEL_H = 1.2 * cm          # vertical space for (a)(b)(c)
BOTTOM_MARGIN = 6 * cm      # caption + breathing room

EXCLUDE_GLACIERS = {103, 104}

# =============================================================================
# METADATA
# =============================================================================
def load_meta():
    return pd.read_csv(CSV_PATH)

def glacier_id(fid, meta):
    return int(meta.loc[meta.feature_ID == fid, "glacier_ID"].iloc[0])

def glacier_name(fid, meta):
    return meta.loc[meta.feature_ID == fid, "glacier_name"].iloc[0]

# =============================================================================
# RASTERISATION
# =============================================================================
def pdf_to_png(pdf):
    png = PNG_CACHE / f"{pdf.stem}.png"

    if png.exists():
        try:
            Image.open(png).verify()
            return png
        except Exception:
            png.unlink(missing_ok=True)

    try:
        img = convert_from_path(pdf, dpi=PNG_DPI, first_page=1, last_page=1)[0]
        img.save(png, "PNG")
        return png
    except Exception as e:
        print(f"[WARNING] Skipping {pdf.name}: {e}")
        return None

# =============================================================================
# PAGE HELPERS
# =============================================================================
def make_page(draw_fn):
    buf = BytesIO()
    c = canvas.Canvas(buf, pagesize=A4)
    draw_fn(c)
    c.showPage()
    c.save()
    buf.seek(0)
    return PdfReader(buf).pages[0]

def draw_page_number(c, n):
    c.setFont("Times-Roman", 10)
    y = PAGE_H - 1.5 * cm
    if n % 2 == 1:
        c.drawRightString(PAGE_W - 2 * cm, y, str(n))
    else:
        c.drawString(2 * cm, y, str(n))

def draw_png_fit_width(c, png, y_top, max_h):
    img = Image.open(png)
    iw, ih = img.size

    scale = PAGE_W / iw
    w = PAGE_W
    h = ih * scale

    if h > max_h:
        scale = max_h / ih
        h = max_h
        w = iw * scale

    c.drawImage(
        str(png),
        (PAGE_W - w) / 2,
        y_top - h,
        w,
        h,
        preserveAspectRatio=True,
        mask="auto",
    )

# =============================================================================
# ASI1
# =============================================================================
def build_ASI1(meta):
    files = []
    for p in ASI1_INPUT.glob("fl_*_picked_lines_chad_v3.pdf"):
        fid = int(re.search(r"fl_(\d+)_", p.name).group(1))
        if fid in EXCLUDE_GLACIERS:
            continue
        files.append((glacier_id(fid, meta), fid, p))

    files.sort(key=lambda x: x[0])

    writer = PdfWriter()
    writer.add_page(PdfReader(ASI_FRONT[1]).pages[0])

    fig, page_no = 11, 2

    for _, fid, pdf in tqdm(files, desc="ASI1"):
        png = pdf_to_png(pdf)
        name = glacier_name(fid, meta)

        page = make_page(lambda c: (
            png and draw_png_fit_width(c, png, PAGE_H - 2 * cm, PAGE_H - BOTTOM_MARGIN - 2 * cm),
            c.setFont("Times-Bold", 11),
            c.drawString(2 * cm, 2.8 * cm, f"Figure S{fig}."),
            c.setFont("Times-Roman", 11),
            c.drawString(6 * cm, 2.8 * cm, f"Generated flowlines at {name} glacier catchment."),
            draw_page_number(c, page_no),
        ))

        writer.add_page(page)
        fig += 1
        page_no += 1

    writer.write(open(OUTFILES[1], "wb"))

# =============================================================================
# ASI2 / ASI3
# =============================================================================
def build_multi(asi, meta, inputs, prefix, caption):
    entries = set()
    for p in inputs["fric1"].glob("*.pdf"):
        m = re.search(r"gl_(\d+)_(\d+)_", p.name)
        fid = int(m.group(1))
        if fid in EXCLUDE_GLACIERS:
            continue
        entries.add((fid, int(m.group(2))))

    entries = sorted(entries, key=lambda x: (glacier_id(x[0], meta), x[1]))

    writer = PdfWriter()
    writer.add_page(PdfReader(ASI_FRONT[asi]).pages[0])

    fig = 74 if asi == 2 else 205
    page_no = 2

    usable_h = PAGE_H - BOTTOM_MARGIN
    panel_h = usable_h / 3
    figure_h = panel_h - LABEL_H

    for fid, fl in tqdm(entries, desc=f"ASI{asi}"):
        name = glacier_name(fid, meta)

        pngs = {}
        for fric in ["fric1", "fric3", "fric2"]:
            pdf = inputs[fric] / f"gl_{fid}_{fl}_{prefix}_{fric}.pdf"
            if pdf.exists():
                pngs[fric] = pdf_to_png(pdf)

        if not pngs:
            continue

        page = make_page(lambda c: (
            # (a)
            c.setFont("Times-Roman", 10),
            c.drawString(2 * cm, PAGE_H - 2 * cm - LABEL_H + 0.2 * cm, "(a)"),
            pngs.get("fric1") and draw_png_fit_width(
                c, pngs["fric1"], PAGE_H - 2 * cm - LABEL_H, figure_h
            ),

            # (b)
            c.drawString(2 * cm, PAGE_H - 2 * cm - panel_h - LABEL_H + 0.2 * cm, "(b)"),
            pngs.get("fric3") and draw_png_fit_width(
                c, pngs["fric3"], PAGE_H - 2 * cm - panel_h - LABEL_H, figure_h
            ),

            # (c)
            c.drawString(2 * cm, PAGE_H - 2 * cm - 2 * panel_h - LABEL_H + 0.2 * cm, "(c)"),
            pngs.get("fric2") and draw_png_fit_width(
                c, pngs["fric2"], PAGE_H - 2 * cm - 2 * panel_h - LABEL_H, figure_h
            ),

            # caption
            c.setFont("Times-Bold", 11),
            c.drawString(2 * cm, 2.8 * cm, f"Figure S{fig}."),
            c.setFont("Times-Roman", 11),
            c.drawString(6 * cm, 2.8 * cm, caption.format(name=name, flowline=fl)),

            draw_page_number(c, page_no),
        ))

        writer.add_page(page)
        fig += 1
        page_no += 1

    writer.write(open(OUTFILES[asi], "wb"))

# =============================================================================
# MAIN
# =============================================================================
if __name__ == "__main__":
    meta = load_meta()

    build_ASI1(meta)

    build_multi(
        2, meta, ASI2_INPUT, "omr",
        "Observed, modelled, and residual velocity anomalies at {name}, flowline {flowline}.",
    )

    build_multi(
        3, meta, ASI3_INPUT, "hm",
        "Similarity index values for {name}, flowline {flowline}.",
    )