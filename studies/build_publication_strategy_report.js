const fs = require("fs");
const path = require("path");
const {
  AlignmentType,
  BorderStyle,
  Document,
  ExternalHyperlink,
  Footer,
  Header,
  HeadingLevel,
  PageBreak,
  PageNumber,
  Paragraph,
  ShadingType,
  Table,
  TableCell,
  TableOfContents,
  TableRow,
  TextRun,
  WidthType,
  Packer,
} = require("docx");

const projectRoot = "E:\\CoWork\\PROJECTS\\PyPSA";
const sourcePath = path.join(
  projectRoot,
  "studies",
  "chapter2_publication_strategy_and_timeline_report.md"
);
const outputPath = path.join(
  projectRoot,
  "studies",
  "chapter2_publication_strategy_and_timeline_report.docx"
);

const NAVY = "17324D";
const TEAL = "247B7B";
const PALE_TEAL = "E8F1EF";
const PALE_BLUE = "EAF0F5";
const PALE_GREY = "F3F5F6";
const MID_GREY = "69747E";
const AMBER = "FFF2CC";
const WHITE = "FFFFFF";

function stripMarkdown(text) {
  return text
    .replace(/\*\*/g, "")
    .replace(/`/g, "")
    .replace(/\\_/g, "_")
    .trim();
}

function inlineRuns(text, options = {}) {
  const runs = [];
  const regex = /(\*\*[^*]+\*\*|`[^`]+`|https?:\/\/[^\s)]+)/g;
  let cursor = 0;
  for (const match of text.matchAll(regex)) {
    if (match.index > cursor) {
      runs.push(
        new TextRun({
          text: text.slice(cursor, match.index),
          size: options.size || 21,
          color: options.color || "20272D",
          italics: options.italics || false,
        })
      );
    }
    const token = match[0];
    if (token.startsWith("**")) {
      runs.push(
        new TextRun({
          text: token.slice(2, -2),
          bold: true,
          size: options.size || 21,
          color: options.color || "20272D",
        })
      );
    } else if (token.startsWith("`")) {
      runs.push(
        new TextRun({
          text: token.slice(1, -1),
          font: "Consolas",
          size: 18,
          color: "34495E",
          shading: { type: ShadingType.CLEAR, fill: PALE_GREY },
        })
      );
    } else {
      runs.push(
        new ExternalHyperlink({
          link: token,
          children: [
            new TextRun({
              text: token,
              style: "Hyperlink",
              size: options.size || 20,
            }),
          ],
        })
      );
    }
    cursor = match.index + token.length;
  }
  if (cursor < text.length) {
    runs.push(
      new TextRun({
        text: text.slice(cursor),
        size: options.size || 21,
        color: options.color || "20272D",
        italics: options.italics || false,
      })
    );
  }
  return runs.length
    ? runs
    : [new TextRun({ text, size: options.size || 21 })];
}

function bodyParagraph(text, options = {}) {
  return new Paragraph({
    children: inlineRuns(text, options),
    spacing: { after: 110, line: 276 },
    alignment: options.alignment || AlignmentType.JUSTIFIED,
    shading: options.shading
      ? { type: ShadingType.CLEAR, fill: options.shading }
      : undefined,
    indent: options.indent,
    border: options.border,
    keepLines: true,
  });
}

function heading(text, level) {
  const map = {
    1: HeadingLevel.HEADING_1,
    2: HeadingLevel.HEADING_2,
    3: HeadingLevel.HEADING_3,
  };
  return new Paragraph({
    text: stripMarkdown(text),
    heading: map[level],
    spacing: { before: level === 1 ? 260 : 190, after: 100 },
    keepNext: true,
  });
}

function parseTable(lines) {
  const rows = lines.map((line) =>
    line
      .trim()
      .replace(/^\|/, "")
      .replace(/\|$/, "")
      .split("|")
      .map((cell) => stripMarkdown(cell))
  );
  if (rows.length > 1 && rows[1].every((cell) => /^:?-{3,}:?$/.test(cell))) {
    rows.splice(1, 1);
  }
  const columnCount = Math.max(...rows.map((row) => row.length));
  const tableWidth = 10100;
  const width = Math.floor(tableWidth / columnCount);
  return new Table({
    width: { size: tableWidth, type: WidthType.DXA },
    columnWidths: Array(columnCount).fill(width),
    alignment: AlignmentType.CENTER,
    rows: rows.map(
      (row, rowIndex) =>
        new TableRow({
          tableHeader: rowIndex === 0,
          cantSplit: true,
          children: Array.from({ length: columnCount }, (_, colIndex) => {
            const value = row[colIndex] || "";
            return new TableCell({
              width: { size: width, type: WidthType.DXA },
              verticalAlign: "center",
              margins: { top: 80, bottom: 80, left: 90, right: 90 },
              shading: {
                type: ShadingType.CLEAR,
                fill:
                  rowIndex === 0
                    ? NAVY
                    : rowIndex % 2 === 0
                      ? PALE_GREY
                      : WHITE,
              },
              borders: {
                top: { style: BorderStyle.SINGLE, size: 2, color: "C7D0D8" },
                bottom: {
                  style: BorderStyle.SINGLE,
                  size: 2,
                  color: "C7D0D8",
                },
                left: { style: BorderStyle.SINGLE, size: 2, color: "C7D0D8" },
                right: { style: BorderStyle.SINGLE, size: 2, color: "C7D0D8" },
              },
              children: [
                new Paragraph({
                  children: inlineRuns(value, {
                    size: columnCount >= 5 ? 15 : 17,
                    color: rowIndex === 0 ? WHITE : "20272D",
                  }).map((run) => {
                    if (rowIndex === 0 && run instanceof TextRun) {
                      return new TextRun({
                        text: run.options?.text || value,
                        bold: true,
                        size: columnCount >= 5 ? 15 : 17,
                        color: WHITE,
                      });
                    }
                    return run;
                  }),
                  spacing: { after: 0, line: 220 },
                }),
              ],
            });
          }),
        })
    ),
  });
}

function buildContent(markdown) {
  const lines = markdown.split(/\r?\n/);
  const children = [];
  let paragraphBuffer = [];
  let index = 0;
  let firstRuleSeen = false;
  let inCode = false;
  let codeLines = [];
  let inEquation = false;
  let equationLines = [];
  let numberedListCounter = 0;
  let activeNumberReference = null;

  function flushParagraph() {
    if (paragraphBuffer.length) {
      const text = paragraphBuffer.join(" ").trim();
      if (text) children.push(bodyParagraph(text));
      paragraphBuffer = [];
    }
  }

  while (index < lines.length) {
    const raw = lines[index];
    const line = raw.trim();
    const isNumberedLine = /^\d+\.\s+/.test(line);
    if (!isNumberedLine) {
      activeNumberReference = null;
    }

    if (line.startsWith("```")) {
      flushParagraph();
      if (inCode) {
        children.push(
          new Paragraph({
            children: [
              new TextRun({
                text: codeLines.join("\n"),
                font: "Consolas",
                size: 17,
                color: "243746",
              }),
            ],
            spacing: { before: 80, after: 120, line: 220 },
            shading: { type: ShadingType.CLEAR, fill: PALE_GREY },
            indent: { left: 260, right: 260 },
          })
        );
        codeLines = [];
        inCode = false;
      } else {
        inCode = true;
      }
      index += 1;
      continue;
    }
    if (inCode) {
      codeLines.push(raw);
      index += 1;
      continue;
    }

    if (line === "\\[") {
      flushParagraph();
      inEquation = true;
      equationLines = [];
      index += 1;
      continue;
    }
    if (line === "\\]" && inEquation) {
      children.push(
        new Paragraph({
          children: [
            new TextRun({
              text: equationLines.join(" "),
              font: "Cambria Math",
              size: 21,
            }),
          ],
          alignment: AlignmentType.CENTER,
          spacing: { before: 100, after: 130 },
        })
      );
      inEquation = false;
      equationLines = [];
      index += 1;
      continue;
    }
    if (inEquation) {
      equationLines.push(line);
      index += 1;
      continue;
    }

    if (!line) {
      flushParagraph();
      index += 1;
      continue;
    }

    if (line === "---") {
      flushParagraph();
      if (!firstRuleSeen) {
        firstRuleSeen = true;
        children.push(new Paragraph({ children: [new PageBreak()] }));
        children.push(
          new Paragraph({
            text: "Contents",
            heading: HeadingLevel.HEADING_1,
            spacing: { after: 180 },
          })
        );
        children.push(
          new TableOfContents("Table of Contents", {
            hyperlink: true,
            headingStyleRange: "1-3",
          })
        );
        children.push(new Paragraph({ children: [new PageBreak()] }));
      }
      index += 1;
      continue;
    }

    if (line.startsWith("|") && index + 1 < lines.length) {
      const next = lines[index + 1].trim();
      if (next.startsWith("|") && /---/.test(next)) {
        flushParagraph();
        const tableLines = [];
        while (index < lines.length && lines[index].trim().startsWith("|")) {
          tableLines.push(lines[index]);
          index += 1;
        }
        children.push(parseTable(tableLines));
        children.push(new Paragraph({ spacing: { after: 80 } }));
        continue;
      }
    }

    if (line.startsWith("# ")) {
      flushParagraph();
      children.push(
        new Paragraph({
          children: [
            new TextRun({
              text: stripMarkdown(line.slice(2)),
              bold: true,
              color: NAVY,
              size: 38,
              font: "Aptos Display",
            }),
          ],
          spacing: { before: 900, after: 220 },
          border: {
            bottom: { style: BorderStyle.SINGLE, size: 12, color: TEAL },
          },
        })
      );
      index += 1;
      continue;
    }
    if (line.startsWith("## ")) {
      flushParagraph();
      const text = line.slice(3);
      if (/^Appendix [A-Z]/.test(text)) {
        children.push(new Paragraph({ children: [new PageBreak()] }));
      }
      children.push(heading(text, 1));
      index += 1;
      continue;
    }
    if (line.startsWith("### ")) {
      flushParagraph();
      children.push(heading(line.slice(4), 2));
      index += 1;
      continue;
    }
    if (line.startsWith("#### ")) {
      flushParagraph();
      children.push(heading(line.slice(5), 3));
      index += 1;
      continue;
    }

    if (/^\*\*[^*]+:\*\*/.test(line) && !firstRuleSeen) {
      flushParagraph();
      children.push(
        bodyParagraph(line, {
          alignment: AlignmentType.LEFT,
          color: MID_GREY,
          size: 19,
        })
      );
      index += 1;
      continue;
    }

    if (/^\d+\.\s+/.test(line)) {
      flushParagraph();
      if (!activeNumberReference) {
        numberedListCounter += 1;
        activeNumberReference = `numbered-list-${numberedListCounter}`;
      }
      children.push(
        new Paragraph({
          children: inlineRuns(line.replace(/^\d+\.\s+/, "")),
          numbering: { reference: activeNumberReference, level: 0 },
          spacing: { after: 70, line: 260 },
        })
      );
      index += 1;
      continue;
    }

    if (/^-\s+/.test(line)) {
      flushParagraph();
      const item = line.replace(/^-\s+/, "");
      const checkbox = /^\[[ xX]\]\s+/.test(item);
      children.push(
        new Paragraph({
          children: checkbox
            ? [
                new TextRun({ text: "☐ ", size: 21, color: TEAL }),
                ...inlineRuns(item.replace(/^\[[ xX]\]\s+/, "")),
              ]
            : inlineRuns(item),
          numbering: checkbox
            ? undefined
            : { reference: "bullet-list", level: 0 },
          indent: checkbox ? { left: 360, hanging: 180 } : undefined,
          spacing: { after: 65, line: 250 },
        })
      );
      index += 1;
      continue;
    }

    if (line.startsWith("> ")) {
      flushParagraph();
      children.push(
        bodyParagraph(line.slice(2), {
          italics: true,
          shading: PALE_TEAL,
          indent: { left: 320, right: 320 },
          border: {
            left: { style: BorderStyle.SINGLE, size: 12, color: TEAL },
          },
        })
      );
      index += 1;
      continue;
    }

    paragraphBuffer.push(line);
    index += 1;
  }

  flushParagraph();
  return children;
}

const markdown = fs.readFileSync(sourcePath, "utf8");
const content = buildContent(markdown);

const doc = new Document({
  creator: "BC-PyPSA research workflow",
  title: "Publication Strategy, Scientific Quality Plan, and Delivery Timeline",
  subject: "Hydropower cascade representation and operational feasibility",
  description:
    "A complete publication-readiness strategy and 20-week delivery plan for Chapter 2.",
  styles: {
    default: {
      document: {
        run: { font: "Aptos", size: 21, color: "20272D" },
        paragraph: { spacing: { after: 100, line: 276 } },
      },
      heading1: {
        run: { font: "Aptos Display", size: 30, bold: true, color: NAVY },
        paragraph: { spacing: { before: 260, after: 110 }, keepNext: true },
      },
      heading2: {
        run: { font: "Aptos Display", size: 25, bold: true, color: TEAL },
        paragraph: { spacing: { before: 210, after: 90 }, keepNext: true },
      },
      heading3: {
        run: { font: "Aptos", size: 22, bold: true, color: NAVY },
        paragraph: { spacing: { before: 160, after: 70 }, keepNext: true },
      },
    },
  },
  numbering: {
    config: [
      {
        reference: "bullet-list",
        levels: [
          {
            level: 0,
            format: "bullet",
            text: "•",
            alignment: AlignmentType.LEFT,
            style: {
              paragraph: { indent: { left: 480, hanging: 240 } },
            },
          },
        ],
      },
      ...Array.from({ length: 60 }, (_, index) => ({
        reference: `numbered-list-${index + 1}`,
        levels: [
          {
            level: 0,
            format: "decimal",
            text: "%1.",
            alignment: AlignmentType.LEFT,
            style: {
              paragraph: { indent: { left: 500, hanging: 260 } },
            },
          },
        ],
      })),
    ],
  },
  sections: [
    {
      properties: {
        page: {
          size: { width: 12240, height: 15840 },
          margin: {
            top: 1080,
            right: 1080,
            bottom: 1000,
            left: 1080,
            header: 400,
            footer: 420,
          },
        },
      },
      headers: {
        default: new Header({
          children: [
            new Paragraph({
              children: [
                new TextRun({
                  text: "BC-PyPSA · Chapter 2 publication strategy",
                  size: 17,
                  color: MID_GREY,
                  italics: true,
                }),
              ],
              border: {
                bottom: { style: BorderStyle.SINGLE, size: 4, color: "D6DCE1" },
              },
              spacing: { after: 50 },
            }),
          ],
        }),
      },
      footers: {
        default: new Footer({
          children: [
            new Paragraph({
              alignment: AlignmentType.RIGHT,
              children: [
                new TextRun({ text: "Page ", size: 17, color: MID_GREY }),
                new TextRun({
                  children: [PageNumber.CURRENT],
                  size: 17,
                  color: MID_GREY,
                }),
              ],
            }),
          ],
        }),
      },
      children: content,
    },
  ],
});

Packer.toBuffer(doc).then((buffer) => {
  fs.writeFileSync(outputPath, buffer);
  console.log(`Wrote ${outputPath}`);
});
