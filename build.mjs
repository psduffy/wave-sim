// build.mjs — inline the ES modules into one self-contained HTML file.
//
// ES modules cannot be loaded from file:// in most browsers, so the standalone
// build wraps each module in a function scope and resolves imports through a
// small registry.  Output: dist/wave-sim.html (open it directly from disk).

import { readFileSync, writeFileSync, mkdirSync } from 'node:fs';

const ORDER = ['physics', 'bathymetry', 'vessel', 'wavefield', 'slowmodel', 'render', 'main'];
const parts = ['const __m = {};'];
for (const name of ORDER) {
  const src = readFileSync(new URL(`./src/${name}.js`, import.meta.url), 'utf8');
  const exports = [...src.matchAll(/^export (?:async )?(?:function|const|let|class) (\w+)/gm)].map((m) => m[1]);
  let body = src
    .replace(/^import \{([^}]*)\} from '\.\/(\w+)\.js';/gm, (_, names, mod) => `const {${names}} = __m.${mod};`)
    .replace(/^import \* as (\w+) from '\.\/(\w+)\.js';/gm, (_, alias, mod) => `const ${alias} = __m.${mod};`)
    .replace(/^export /gm, '');
  parts.push(`__m.${name} = (function () {\n${body}\nreturn { ${exports.join(', ')} };\n})();`);
}
const bundle = parts.join('\n\n');

const page = readFileSync(new URL('./index.html', import.meta.url), 'utf8');
const tag = '<script type="module" src="./src/main.js"></script>';
if (!page.includes(tag)) throw new Error('module script tag not found in index.html');
const standalone = page.replace(tag, `<script>\n${bundle}\n</script>`);
mkdirSync(new URL('./dist/', import.meta.url), { recursive: true });
writeFileSync(new URL('./dist/wave-sim.html', import.meta.url), standalone);

// Variant without the document skeleton, for hosts that wrap the page themselves.
const fragment = standalone
  .replace(/^<!doctype html>\s*/i, '')
  .replace(/<html[^>]*>\s*/i, '')
  .replace(/<\/html>\s*$/i, '')
  .replace(/<head>\s*/i, '')
  .replace(/<\/head>\s*/i, '')
  .replace(/<body>\s*/i, '')
  .replace(/<\/body>\s*/i, '')
  .replace(/<meta charset="utf-8">\s*/i, '')
  .replace(/<meta name="viewport"[^>]*>\s*/i, '');
writeFileSync(new URL('./dist/wave-sim.fragment.html', import.meta.url), fragment);
console.log(`built dist/wave-sim.html (${(standalone.length / 1024).toFixed(0)} kB) and dist/wave-sim.fragment.html`);
