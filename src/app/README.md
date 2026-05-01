# App icon resources

`sketcher.icns` (macOS), `sketcher.ico` (Windows), and `sketcher.rc` (Windows
resource script) are baked into the standalone `Sketcher` executable. The
`.icns` and `.ico` files are generated from the master SVG at
`src/schrodinger/sketcher/icons/sketcher-logo.svg` and should be regenerated
whenever that SVG changes.

## Prerequisites

ImageMagick and librsvg. On macOS:

```sh
brew install imagemagick librsvg
```

`iconutil` ships with macOS.

> **Important:** rasterize with `rsvg-convert` (from `librsvg`), not directly
> with ImageMagick. The logo uses linear gradients that ImageMagick's built-in
> MSVG renderer mishandles, producing wrong colors. Going through
> `rsvg-convert` for the PNGs and only using ImageMagick for the final `.ico`
> container assembly avoids the problem.

## Regenerate `sketcher.ico` (Windows)

```sh
SRC=src/schrodinger/sketcher/icons/sketcher-logo.svg
TMP=$(mktemp -d)
for size in 16 24 32 48 64 128 256; do
  rsvg-convert -w $size -h $size "$SRC" -o "$TMP/icon_${size}.png"
done
magick "$TMP"/icon_{16,24,32,48,64,128,256}.png src/app/sketcher.ico
```

## Regenerate `sketcher.icns` (macOS)

```sh
SRC=src/schrodinger/sketcher/icons/sketcher-logo.svg
ICONSET=$(mktemp -d)/sketcher.iconset
mkdir -p "$ICONSET"
for spec in 16:icon_16x16.png 32:icon_16x16@2x.png \
            32:icon_32x32.png 64:icon_32x32@2x.png \
            128:icon_128x128.png 256:icon_128x128@2x.png \
            256:icon_256x256.png 512:icon_256x256@2x.png \
            512:icon_512x512.png 1024:icon_512x512@2x.png; do
  size=${spec%%:*}
  name=${spec##*:}
  rsvg-convert -w $size -h $size "$SRC" -o "$ICONSET/$name"
done
iconutil -c icns "$ICONSET" -o src/app/sketcher.icns
```
