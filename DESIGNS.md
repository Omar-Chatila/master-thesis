# Thesis design variants

All variants compile the same `main.tex`; only typography, colour, cover artwork, and chapter-number treatment change.

1. **Habitat Grid** (`habitat-grid`) — the original balanced crossover design: ecological palette, analytical grid, and trajectory curve.
2. **Field Notes** (`field-notes`) — warmer and more book-like, with an olive field-journal palette and topographic contours.
3. **Technical Atlas** (`technical-atlas`) — crisp blue/teal geometry emphasizing GIS and algorithm engineering.
4. **Quiet Academic** (`quiet-academic`) — restrained and formal, with minimal decoration and serif chapter titles.
5. **Image Watermark** (`image-watermark`) — a genuinely different centered, classical layout with a replaceable translucent PNG behind the title-page information and centered chapter openers.
6. **Tengwar Book** (`tengwar-book`) — a two-stage dark outer cover and sparse inner title page, EB Garamond typography, centered word-number chapter openers, and book-style running folios.

## Replacing the title-page PNG

`Image Watermark` uses `figures/title-watermark-placeholder.png`. Replace that file with your own PNG, or change `\ThesisTitleImage` in `config/thesis-meta.tex`. Adjust `\ThesisTitleImageOpacity` there; `0.10` is subtle and `0.20` is comparatively strong. The image is constrained to the central title-page area and keeps its aspect ratio.

## Tengwar attribution

`Tengwar Book` is an independent adaptation for this repository based on the layout principles of Níckolas de Aguiar Alves's [Tengwar Thesis Template](https://www.overleaf.com/latex/templates/tengwar-thesis-template/xpkdstnvfsgg), published under CC BY 4.0. It retains this project's KOMA-Script structure and thesis content; it does not copy or require the original custom class. The included EB Garamond font files are distributed under the SIL Open Font License in `assets/fonts/OFL.txt`.

Build all complete variants with:

```bash
make designs
```

The PDFs are written to `build/designs/`. All variants compile with pdfLaTeX; `tengwar-book` additionally requires the standard `ebgaramond` TeX package. To select one as the default, change `\ThesisDesign` in `config/thesis-meta.tex` and run `make`.
