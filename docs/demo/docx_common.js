// Общие стили и помощники для документов Word (docx-js).
const fs = require("fs");
const path = require("path");
const d = require("docx");

const FONT = "Arial";
const W = 9638; // ширина текста A4 при полях 1134 DXA (2 см): 11906 − 2·1134

const styles = {
  default: { document: { run: { font: FONT, size: 21 }, paragraph: { spacing: { after: 100, line: 276 } } } },
  paragraphStyles: [
    { id: "Title", name: "Title", basedOn: "Normal", next: "Normal", quickFormat: true,
      run: { font: FONT, size: 44, bold: true, color: "1F3A5F" }, paragraph: { spacing: { before: 0, after: 200 }, keepNext: true } },
    { id: "Heading1", name: "Heading 1", basedOn: "Normal", next: "Normal", quickFormat: true,
      run: { font: FONT, size: 30, bold: true, color: "1F3A5F" }, paragraph: { spacing: { before: 360, after: 140 }, keepNext: true, outlineLevel: 0 } },
    { id: "Heading2", name: "Heading 2", basedOn: "Normal", next: "Normal", quickFormat: true,
      run: { font: FONT, size: 25, bold: true, color: "2A4D7A" }, paragraph: { spacing: { before: 240, after: 100 }, keepNext: true, outlineLevel: 1 } },
    { id: "Heading3", name: "Heading 3", basedOn: "Normal", next: "Normal", quickFormat: true,
      run: { font: FONT, size: 22, bold: true, color: "2A4D7A" }, paragraph: { spacing: { before: 180, after: 80 }, keepNext: true, outlineLevel: 2 } },
    { id: "Caption", name: "Caption", basedOn: "Normal", next: "Normal", quickFormat: true,
      run: { font: FONT, size: 18, italics: true, color: "444444" }, paragraph: { spacing: { before: 60, after: 200 } } },
    { id: "Note", name: "Note", basedOn: "Normal", next: "Normal",
      run: { font: FONT, size: 19, color: "333333" }, paragraph: { spacing: { before: 80, after: 120 }, indent: { left: 360 } } },
  ],
};

const numbering = {
  config: [
    { reference: "bullets", levels: [{ level: 0, format: d.LevelFormat.BULLET, text: "•", alignment: d.AlignmentType.LEFT,
        style: { paragraph: { indent: { left: 540, hanging: 270 } } } },
      { level: 1, format: d.LevelFormat.BULLET, text: "–", alignment: d.AlignmentType.LEFT,
        style: { paragraph: { indent: { left: 1000, hanging: 270 } } } }] },
    ...Array.from({ length: 40 }, (_, i) => ({ reference: "steps" + i, levels: [{ level: 0, format: d.LevelFormat.DECIMAL,
      text: "%1.", alignment: d.AlignmentType.LEFT, style: { paragraph: { indent: { left: 540, hanging: 360 } } } }] })),
  ],
};

// Текст с простой разметкой: **жирный**, *курсив*, `моноширинный`
function runs(text, base = {}) {
  const out = [];
  const re = /(\*\*[^*]+\*\*|\*[^*]+\*|`[^`]+`)/g;
  let last = 0, m;
  while ((m = re.exec(text))) {
    if (m.index > last) out.push(new d.TextRun({ text: text.slice(last, m.index), ...base }));
    const t = m[0];
    if (t.startsWith("**")) out.push(new d.TextRun({ text: t.slice(2, -2), bold: true, ...base }));
    else if (t.startsWith("`")) out.push(new d.TextRun({ text: t.slice(1, -1), font: "Consolas", size: 19, ...base }));
    else out.push(new d.TextRun({ text: t.slice(1, -1), italics: true, ...base }));
    last = m.index + t.length;
  }
  if (last < text.length) out.push(new d.TextRun({ text: text.slice(last), ...base }));
  return out;
}

const P = (text, opt = {}) => new d.Paragraph({ children: runs(text), ...opt });
const H1 = (t) => new d.Paragraph({ heading: d.HeadingLevel.HEADING_1, children: [new d.TextRun(t)] });
const H2 = (t) => new d.Paragraph({ heading: d.HeadingLevel.HEADING_2, children: [new d.TextRun(t)] });
const H3 = (t) => new d.Paragraph({ heading: d.HeadingLevel.HEADING_3, children: [new d.TextRun(t)] });
const B = (t, level = 0) => new d.Paragraph({ numbering: { reference: "bullets", level }, children: runs(t) });
let stepList = 0;
function Steps(items) {
  const ref = "steps" + (stepList++);
  return items.map((t) => new d.Paragraph({ numbering: { reference: ref, level: 0 }, children: runs(t), spacing: { after: 80 } }));
}
const Note = (label, text, color = "1F3A5F") =>
  new d.Paragraph({ style: "Note", border: { left: { style: d.BorderStyle.SINGLE, size: 12, color, space: 8 } },
    children: [new d.TextRun({ text: label + "  ", bold: true, color }), ...runs(text)] });
const Code = (t) => new d.Paragraph({ shading: { type: d.ShadingType.CLEAR, fill: "F3F3F0", color: "auto" },
  spacing: { after: 120 }, indent: { left: 240 }, children: [new d.TextRun({ text: t, font: "Consolas", size: 18 })] });

function Fig(file, widthPx, caption) {
  const buf = fs.readFileSync(file);
  // размеры PNG из заголовка IHDR
  const w = buf.readUInt32BE(16), h = buf.readUInt32BE(20);
  const height = Math.round(widthPx * h / w);
  return [
    new d.Paragraph({ alignment: d.AlignmentType.CENTER, keepNext: true, spacing: { before: 120, after: 40 },
      children: [new d.ImageRun({ type: "png", data: buf, transformation: { width: widthPx, height },
        altText: { title: path.basename(file), description: caption, name: path.basename(file) } })] }),
    new d.Paragraph({ style: "Caption", alignment: d.AlignmentType.CENTER, children: runs(caption) }),
  ];
}

const border = { style: d.BorderStyle.SINGLE, size: 4, color: "BFBFBF" };
const borders = { top: border, bottom: border, left: border, right: border };
function Table(header, rows, widths) {
  const total = widths.reduce((a, b) => a + b, 0);
  const cell = (t, i, head) => new d.TableCell({
    borders, width: { size: widths[i], type: d.WidthType.DXA },
    shading: head ? { fill: "DCE6F2", type: d.ShadingType.CLEAR, color: "auto" } : undefined,
    margins: { top: 60, bottom: 60, left: 100, right: 100 },
    children: String(t).split("\n").map((line) => new d.Paragraph({ spacing: { after: 0 },
      children: runs(line, head ? { bold: true, size: 19 } : { size: 19 }) })),
  });
  return new d.Table({
    width: { size: total, type: d.WidthType.DXA }, columnWidths: widths,
    rows: [new d.TableRow({ tableHeader: true, children: header.map((h, i) => cell(h, i, true)) }),
      ...rows.map((r) => new d.TableRow({ children: r.map((c, i) => cell(c, i, false)) }))],
  });
}
const Spacer = () => new d.Paragraph({ spacing: { after: 60 }, children: [] });

function footer(text) {
  return { default: new d.Footer({ children: [new d.Paragraph({ alignment: d.AlignmentType.CENTER,
    children: [new d.TextRun({ text: text + "   ·   ", size: 16, color: "777777" }),
      new d.TextRun({ children: [d.PageNumber.CURRENT], size: 16, color: "777777" })] })] }) };
}
const page = { size: { width: 11906, height: 16838 }, margin: { top: 1134, bottom: 1134, left: 1134, right: 1134 } };

async function write(doc, file) {
  const buf = await d.Packer.toBuffer(doc);
  fs.writeFileSync(file, buf);
  console.log("written", file, buf.length);
}

module.exports = { d, styles, numbering, runs, P, H1, H2, H3, B, Steps, Note, Code, Fig, Table, Spacer, footer, page, write, W };
