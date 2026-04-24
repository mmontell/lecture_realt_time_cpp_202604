# Scientific Lecture Template

A minimalist HTML slide template for computational science × high-energy
physics talks. White / blue / black palette, abstract dots-and-lines motifs,
KaTeX math, syntax-highlighted code, and a built-in Tweaks panel.

**Footer credit:** Marco Montella — Ohio State University — Guest Lecture at
University of Manchester, 2026.

---

## Files in this project

| File | What it is |
|---|---|
| `Scientific Lecture Template.html` | **The deck itself.** Open this. |
| `deck.css` | All styling (type, colors, slide layouts, dark mode). |
| `deck-stage.js` | Slide navigation, scaling, keyboard shortcuts, speaker notes bridge. |
| `motifs.js` | Procedural SVG generators for the detector / graph / lattice motifs. |

Everything is plain text. No build step. No bundler. Refresh the browser to
see changes.

---

## How to present

1. Open `Scientific Lecture Template.html` in any modern browser (Chrome,
   Firefox, Safari — doesn't matter).
2. Navigate with:
   - **← / →** — prev / next slide
   - **Home / End** — first / last slide
   - **F** — fullscreen
   - **P** or **Ctrl/Cmd + P** — print / save as PDF (one page per slide)
   - **N** — toggle speaker-notes window
3. Your position is auto-saved to `localStorage`, so refreshing doesn't lose
   your place.

---

## How to edit content

Open the `.html` file in any text editor (VS Code, Sublime, even TextEdit).
Each slide is a `<section>` inside `<deck-stage>`. Find the slide you want
and edit its contents directly. Save → refresh the browser.

### Titles and bullets

Just plain HTML:

```html
<h1>My Slide Title</h1>
<ul>
  <li>First point</li>
  <li>Second point</li>
</ul>
```

### Math (KaTeX)

Inline: `\( E = mc^2 \)`
Display: `\[ \mathcal{L} = -\tfrac{1}{4} F_{\mu\nu} F^{\mu\nu} \]`

KaTeX re-renders automatically whenever you advance a slide.

### Code snippets

Wrap in `<pre><code class="language-cpp">…</code></pre>`. Highlight.js
picks up the language class. Supported out of the box: `cpp`, `python`,
`bash`, `json`, plus any language highlight.js supports.

### Speaker notes

There's a single JSON array near the top of the file:

```html
<script type="application/json" id="speaker-notes">
[
  "Notes for slide 1…",
  "Notes for slide 2…"
]
</script>
```

One entry per slide, in order. Press **N** during the talk to pop out the
notes window on a second screen.

### Adding a new slide

Copy any existing `<section class="slide …">…</section>` block inside
`<deck-stage>` and paste it where you want it. Add a matching speaker-notes
string to the JSON array. That's it.

---

## The Tweaks panel

Click the **Tweaks** toggle in the top toolbar (when previewing inside this
editor) to change accent color, font pairing, density, motif on/off, and
dark mode live. Changes persist to the file's `TWEAK_DEFAULTS` block so
they survive reloads.

Outside the editor (plain browser), the panel is hidden — the defaults
already baked into the file are what you'll see.

---

## Exporting for other formats

| You want… | Do this |
|---|---|
| **PDF for sharing** | Ctrl/Cmd + P in the browser → Save as PDF |
| **Present on someone else's laptop** | Send them the whole project folder, or the PDF |
| **Editable PowerPoint** | Ask me to export as PPTX (editable). Caveat: math becomes images, code loses syntax colors, motifs become static pictures |
| **Pixel-perfect PowerPoint** | Ask me to export as PPTX (screenshots). Looks identical, but no text is editable |
| **Self-contained single HTML file** | Ask me to bundle it — one file with all assets inlined, works offline |
| **Google Slides** | No direct path. Export to PPTX first, then upload to Google Slides (same caveats apply) |

**Recommendation:** keep editing the HTML, export a PDF for each talk.
You get the math, code highlighting, motifs, and speaker notes exactly as
designed, and the PDF is what you actually project.

---

## What PPTX / Google Slides can't preserve

Being upfront about the trade-offs if you do convert:

- **LaTeX math** → rasterized to images. Still looks fine, not editable as math.
- **Syntax-highlighted code** → plain monospace text, no colors.
- **Procedural SVG motifs** (detector rings, graph, lattice) → flat images.
  Regenerating with a new random seed means re-exporting.
- **Live Tweaks / dark mode** → frozen to whatever was set at export time.
- **Speaker notes** → preserved (they copy into PowerPoint's notes field).

If any of those matter to you, stay in HTML.

---

## Troubleshooting

- **Math not rendering?** Check the browser console. KaTeX loads from a CDN;
  you need internet on first open.
- **Fonts look wrong?** Same — Google Fonts loads from a CDN.
- **Offline presenting?** Ask me to produce a self-contained bundle; it
  inlines fonts, KaTeX, and highlight.js so nothing depends on network.
- **Slide counter stuck?** Refresh the page; `localStorage` remembered an
  old slide index from a previous version.
