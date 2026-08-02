"""Build the Chapter 2 journal manuscript DOCX from the reproducible Markdown source.

The source proposal is read only to retain its four tables and four embedded
working figures. The generated document is written to studies/build first; it
is copied over the proposal path only after rendering and validation.
"""

from __future__ import annotations

import argparse
import re
from copy import deepcopy
from io import BytesIO
from pathlib import Path

from docx import Document
from docx.enum.section import WD_SECTION
from docx.enum.style import WD_STYLE_TYPE
from docx.enum.table import WD_CELL_VERTICAL_ALIGNMENT, WD_TABLE_ALIGNMENT
from docx.enum.text import WD_ALIGN_PARAGRAPH, WD_BREAK, WD_LINE_SPACING
from docx.oxml import OxmlElement
from docx.oxml.ns import qn
from docx.shared import Inches, Pt, RGBColor
from PIL import Image


PROJECT_ROOT = Path(r"E:\CoWork\PROJECTS\PyPSA")
DEFAULT_SOURCE = PROJECT_ROOT / "studies" / "chapter2_journal_manuscript_source.md"
DEFAULT_PROPOSAL = Path(
    r"E:\CoWork\PROJECTS\QE\CH_02\2026 07 01 - Chapter2_Cascade_Draft_EL_TN_EL.docx"
)
DEFAULT_OUTPUT = PROJECT_ROOT / "studies" / "build" / "chapter2_journal_manuscript.docx"

NAVY = "16324F"
TEAL = "247B7B"
LIGHT_TEAL = "E8F1EF"
PLACEHOLDER = "FFF2CC"
LIGHT_GREY = "EEF1F3"
DARK_GREY = RGBColor(70, 78, 86)


def set_cell_shading(cell, fill: str) -> None:
    tc_pr = cell._tc.get_or_add_tcPr()
    shd = tc_pr.find(qn("w:shd"))
    if shd is None:
        shd = OxmlElement("w:shd")
        tc_pr.append(shd)
    shd.set(qn("w:fill"), fill)


def set_paragraph_shading(paragraph, fill: str) -> None:
    p_pr = paragraph._p.get_or_add_pPr()
    shd = p_pr.find(qn("w:shd"))
    if shd is None:
        shd = OxmlElement("w:shd")
        p_pr.append(shd)
    shd.set(qn("w:fill"), fill)


def set_cell_margins(cell, top=80, start=90, bottom=80, end=90) -> None:
    tc = cell._tc
    tc_pr = tc.get_or_add_tcPr()
    tc_mar = tc_pr.first_child_found_in("w:tcMar")
    if tc_mar is None:
        tc_mar = OxmlElement("w:tcMar")
        tc_pr.append(tc_mar)
    for margin, value in (("top", top), ("start", start), ("bottom", bottom), ("end", end)):
        node = tc_mar.find(qn(f"w:{margin}"))
        if node is None:
            node = OxmlElement(f"w:{margin}")
            tc_mar.append(node)
        node.set(qn("w:w"), str(value))
        node.set(qn("w:type"), "dxa")


def add_page_number(paragraph) -> None:
    paragraph.alignment = WD_ALIGN_PARAGRAPH.RIGHT
    run = paragraph.add_run("Page ")
    run.font.size = Pt(9)
    begin = OxmlElement("w:fldChar")
    begin.set(qn("w:fldCharType"), "begin")
    instruction = OxmlElement("w:instrText")
    instruction.set(qn("xml:space"), "preserve")
    instruction.text = "PAGE"
    separate = OxmlElement("w:fldChar")
    separate.set(qn("w:fldCharType"), "separate")
    text = OxmlElement("w:t")
    text.text = "1"
    end = OxmlElement("w:fldChar")
    end.set(qn("w:fldCharType"), "end")
    run._r.extend([begin, instruction, separate, text, end])


def configure_document(doc: Document) -> None:
    section = doc.sections[0]
    section.page_width = Inches(8.5)
    section.page_height = Inches(11)
    section.top_margin = Inches(0.82)
    section.bottom_margin = Inches(0.78)
    section.left_margin = Inches(0.9)
    section.right_margin = Inches(0.9)
    section.header_distance = Inches(0.32)
    section.footer_distance = Inches(0.32)

    normal = doc.styles["Normal"]
    normal.font.name = "Aptos"
    normal.font.size = Pt(10.5)
    normal.font.color.rgb = RGBColor(28, 34, 40)
    normal.paragraph_format.space_after = Pt(5.5)
    normal.paragraph_format.line_spacing_rule = WD_LINE_SPACING.SINGLE
    normal.paragraph_format.line_spacing = 1.08

    title = doc.styles["Title"]
    title.font.name = "Aptos Display"
    title.font.size = Pt(20)
    title.font.bold = True
    title.font.color.rgb = RGBColor.from_string(NAVY)
    title.paragraph_format.space_after = Pt(10)

    for style_name, size, before, after in (
        ("Heading 1", 15, 12, 6),
        ("Heading 2", 12.5, 10, 4),
        ("Heading 3", 11, 8, 3),
    ):
        style = doc.styles[style_name]
        style.font.name = "Aptos Display"
        style.font.size = Pt(size)
        style.font.bold = True
        style.font.color.rgb = RGBColor.from_string(NAVY if style_name != "Heading 3" else TEAL)
        style.paragraph_format.space_before = Pt(before)
        style.paragraph_format.space_after = Pt(after)
        style.paragraph_format.keep_with_next = True

    if "Placeholder" not in [s.name for s in doc.styles]:
        style = doc.styles.add_style("Placeholder", WD_STYLE_TYPE.PARAGRAPH)
        style.base_style = doc.styles["Normal"]
        style.font.name = "Aptos"
        style.font.size = Pt(9.5)
        style.font.bold = False
        style.font.color.rgb = RGBColor(90, 67, 12)
        style.paragraph_format.left_indent = Inches(0.15)
        style.paragraph_format.right_indent = Inches(0.15)
        style.paragraph_format.space_before = Pt(4)
        style.paragraph_format.space_after = Pt(6)

    if "Equation" not in [s.name for s in doc.styles]:
        style = doc.styles.add_style("Equation", WD_STYLE_TYPE.PARAGRAPH)
        style.base_style = doc.styles["Normal"]
        style.font.name = "Cambria Math"
        style.font.size = Pt(10)
        style.paragraph_format.left_indent = Inches(0.35)
        style.paragraph_format.right_indent = Inches(0.35)
        style.paragraph_format.space_before = Pt(4)
        style.paragraph_format.space_after = Pt(6)

    if "Manuscript Note" not in [s.name for s in doc.styles]:
        style = doc.styles.add_style("Manuscript Note", WD_STYLE_TYPE.PARAGRAPH)
        style.base_style = doc.styles["Normal"]
        style.font.name = "Aptos"
        style.font.size = Pt(9)
        style.font.italic = True
        style.font.color.rgb = DARK_GREY

    header = section.header.paragraphs[0]
    header.text = "Hydropower cascade aggregation and operational feasibility"
    header.style = doc.styles["Manuscript Note"]
    header.alignment = WD_ALIGN_PARAGRAPH.LEFT
    header.paragraph_format.space_after = Pt(0)
    footer = section.footer.paragraphs[0]
    add_page_number(footer)


def rich_text(paragraph, text: str) -> None:
    """Add lightweight Markdown emphasis without leaving markup in the DOCX."""
    text = text.replace("**", "\x00")
    pieces = text.split("\x00")
    for idx, piece in enumerate(pieces):
        if not piece:
            continue
        # Preserve inline code as plain text in a monospace font.
        code_parts = re.split(r"(`[^`]+`)", piece)
        for cp in code_parts:
            if not cp:
                continue
            run = paragraph.add_run(cp[1:-1] if cp.startswith("`") and cp.endswith("`") else cp)
            run.bold = idx % 2 == 1
            if cp.startswith("`") and cp.endswith("`"):
                run.font.name = "Consolas"
                run.font.size = Pt(9)


def extract_source_material(proposal_path: Path):
    source = Document(proposal_path)
    tables = []
    for table in source.tables:
        rows = []
        for row in table.rows:
            rows.append([cell.text.strip() for cell in row.cells])
        tables.append(rows)

    images = []
    seen = set()
    for shape in source.inline_shapes:
        rid = shape._inline.graphic.graphicData.pic.blipFill.blip.embed
        part = source.part.related_parts[rid]
        key = str(part.partname)
        if key not in seen:
            images.append((key, part.blob))
            seen.add(key)
    return tables, images


def add_source_table(doc: Document, rows: list[list[str]], table_no: int) -> None:
    if not rows:
        return
    max_cols = max(len(r) for r in rows)
    table = doc.add_table(rows=len(rows), cols=max_cols)
    table.alignment = WD_TABLE_ALIGNMENT.CENTER
    table.style = "Table Grid"
    for r_idx, values in enumerate(rows):
        for c_idx in range(max_cols):
            cell = table.cell(r_idx, c_idx)
            cell.vertical_alignment = WD_CELL_VERTICAL_ALIGNMENT.CENTER
            set_cell_margins(cell)
            value = values[c_idx] if c_idx < len(values) else ""
            cell.text = value
            for p in cell.paragraphs:
                p.paragraph_format.space_after = Pt(1.5)
                p.paragraph_format.line_spacing = 1.0
                for run in p.runs:
                    run.font.name = "Aptos"
                    run.font.size = Pt(7.2 if max_cols >= 7 else 8)
                    if r_idx == 0:
                        run.bold = True
                        run.font.color.rgb = RGBColor(255, 255, 255)
            if r_idx == 0:
                set_cell_shading(cell, NAVY)
            elif r_idx % 2 == 0:
                set_cell_shading(cell, LIGHT_GREY)
    # Repeat header row in Word.
    tr_pr = table.rows[0]._tr.get_or_add_trPr()
    tbl_header = OxmlElement("w:tblHeader")
    tbl_header.set(qn("w:val"), "true")
    tr_pr.append(tbl_header)

    p = doc.add_paragraph(style="Manuscript Note")
    p.alignment = WD_ALIGN_PARAGRAPH.CENTER
    p.add_run(
        f"Working Table {table_no} retained from the proposal. "
        "Verify all entries, units, and citations before submission."
    )


def add_working_figure(doc: Document, image_blob: bytes, figure_no: int) -> None:
    p = doc.add_paragraph()
    p.alignment = WD_ALIGN_PARAGRAPH.CENTER
    run = p.add_run()
    try:
        run.add_picture(BytesIO(image_blob), width=Inches(6.55))
    except Exception:
        # Convert proposal vector formats (notably EMF/WMF) to a raster image
        # because python-docx does not recognize their headers.
        source_image = Image.open(BytesIO(image_blob))
        converted = BytesIO()
        source_image.convert("RGB").save(converted, format="PNG", dpi=(220, 220))
        converted.seek(0)
        run.add_picture(converted, width=Inches(6.55))
    caption = doc.add_paragraph(style="Manuscript Note")
    caption.alignment = WD_ALIGN_PARAGRAPH.CENTER
    caption.add_run(
        f"Working Figure {figure_no} retained from the proposal. "
        "Replace or regenerate from versioned analysis outputs before submission."
    )


def flush_buffer(doc: Document, lines: list[str], in_references: bool) -> None:
    if not lines:
        return
    text = " ".join(line.strip() for line in lines).strip()
    if not text:
        lines.clear()
        return
    placeholder = bool(
        re.match(
            r"^\*\*\[(RESULT|MODEL|DATA|IMPLEMENTATION|TABLE|FIGURE|SCENARIO|"
            r"VALIDATION|METRIC|DISCUSSION|CONCLUSION|AVAILABILITY|AUTHOR|"
            r"ACKNOWLEDGEMENT|REFERENCE) PLACEHOLDER",
            text,
        )
        or re.match(r"^\*\*\[(TABLE|FIGURE) \d+ ABOUT HERE", text)
    )
    p = doc.add_paragraph(style="Placeholder" if placeholder else "Normal")
    if placeholder:
        text = text.replace("**", "")
        set_paragraph_shading(p, PLACEHOLDER)
        p.paragraph_format.keep_together = True
    rich_text(p, text)
    if in_references and not placeholder:
        p.paragraph_format.left_indent = Inches(0.28)
        p.paragraph_format.first_line_indent = Inches(-0.28)
        p.paragraph_format.space_after = Pt(3)
        for run in p.runs:
            run.font.size = Pt(9)
    lines.clear()


def build(source_md: Path, proposal_docx: Path, output_docx: Path) -> None:
    tables, images = extract_source_material(proposal_docx)
    doc = Document()
    configure_document(doc)

    lines = source_md.read_text(encoding="utf-8").splitlines()
    buffer: list[str] = []
    in_equation = False
    equation_lines: list[str] = []
    in_references = False
    title_written = False

    table_map = {1: 0, 2: 1, 3: 3, 4: 2}
    figure_map = {1: 0, 2: 1, 6: 2, 7: 3}

    for raw in lines:
        stripped = raw.strip()

        if stripped == r"\[":
            flush_buffer(doc, buffer, in_references)
            in_equation = True
            equation_lines = []
            continue
        if stripped == r"\]" and in_equation:
            p = doc.add_paragraph(style="Equation")
            p.alignment = WD_ALIGN_PARAGRAPH.CENTER
            p.add_run(" ".join(equation_lines))
            in_equation = False
            equation_lines = []
            continue
        if in_equation:
            equation_lines.append(stripped)
            continue

        if not stripped:
            flush_buffer(doc, buffer, in_references)
            continue

        if stripped.startswith("# "):
            flush_buffer(doc, buffer, in_references)
            p = doc.add_paragraph(style="Title")
            p.alignment = WD_ALIGN_PARAGRAPH.LEFT
            p.add_run(stripped[2:].strip())
            title_written = True
            note = doc.add_paragraph(style="Manuscript Note")
            note.add_run(
                "Working journal manuscript. Amber callouts mark required analysis, "
                "validation, metadata, or editorial completion."
            )
            set_paragraph_shading(note, LIGHT_TEAL)
            continue

        heading_match = re.match(r"^(#{2,4})\s+(.*)$", stripped)
        if heading_match:
            flush_buffer(doc, buffer, in_references)
            level = len(heading_match.group(1)) - 1
            heading_text = heading_match.group(2)
            in_references = heading_text == "References"
            p = doc.add_paragraph(style=f"Heading {min(level, 3)}")
            rich_text(p, heading_text)
            continue

        if stripped.startswith("- "):
            flush_buffer(doc, buffer, in_references)
            p = doc.add_paragraph(style="List Bullet")
            rich_text(p, stripped[2:])
            continue

        table_marker = re.search(r"\[TABLE\s+(\d+)\s+ABOUT HERE", stripped)
        figure_marker = re.search(r"\[FIGURE\s+(\d+)\s+ABOUT HERE", stripped)
        if table_marker:
            flush_buffer(doc, buffer, in_references)
            p = doc.add_paragraph(style="Placeholder")
            set_paragraph_shading(p, PLACEHOLDER)
            rich_text(p, stripped.replace("**", ""))
            number = int(table_marker.group(1))
            source_index = table_map.get(number)
            if source_index is not None and source_index < len(tables):
                add_source_table(doc, tables[source_index], number)
            continue
        if figure_marker:
            flush_buffer(doc, buffer, in_references)
            p = doc.add_paragraph(style="Placeholder")
            set_paragraph_shading(p, PLACEHOLDER)
            rich_text(p, stripped.replace("**", ""))
            number = int(figure_marker.group(1))
            source_index = figure_map.get(number)
            if source_index is not None and source_index < len(images):
                add_working_figure(doc, images[source_index][1], number)
            continue

        buffer.append(stripped)

    flush_buffer(doc, buffer, in_references)

    # Document properties and a compact build note.
    props = doc.core_properties
    props.title = "Quantifying Operational Bias from Aggregated Hydropower Reservoirs"
    props.subject = "Journal manuscript on hydropower cascade representation in BC-PyPSA"
    props.keywords = "PyPSA, hydropower, cascade, reservoir aggregation, British Columbia"
    props.comments = (
        "Generated from studies/chapter2_journal_manuscript_source.md; "
        "working figures and tables retained from the source proposal."
    )
    output_docx.parent.mkdir(parents=True, exist_ok=True)
    doc.save(output_docx)
    print(f"Wrote {output_docx}")
    print(f"Retained {len(tables)} proposal tables and {len(images)} proposal figures")


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--source", type=Path, default=DEFAULT_SOURCE)
    parser.add_argument("--proposal", type=Path, default=DEFAULT_PROPOSAL)
    parser.add_argument("--output", type=Path, default=DEFAULT_OUTPUT)
    args = parser.parse_args()
    build(args.source, args.proposal, args.output)


if __name__ == "__main__":
    main()
