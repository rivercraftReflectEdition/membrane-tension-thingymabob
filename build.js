'use strict';
/* Build: inline src/ into a single self-contained index.html (GitHub Pages
 * serves the built file; no external requests at runtime).
 *   node build.js
 * The physics core is wrapped in @physics-core sentinels so the tests can
 * verify that the shipped page contains exactly src/physics-core.js. */
const fs = require('fs');
const path = require('path');
const read = f => fs.readFileSync(path.join(__dirname, f), 'utf8');

const html = read('src/template.html')
  .replace('/*__FONTS__*/', () => read('src/fonts.css').trim())
  .replace('/*__STYLES__*/', () => read('src/styles.css').trim())
  .replace('/*__CORE__*/', () =>
    '/* @physics-core-begin */\n' + read('src/physics-core.js').trim() + '\n/* @physics-core-end */')
  .replace('/*__APP__*/', () => read('src/app.js').trim());

if (html.includes('/*__')) {
  console.error('build FAILED: unreplaced placeholder remains');
  process.exit(1);
}
fs.writeFileSync(path.join(__dirname, 'index.html'), html);
console.log('built index.html (%d KB)', Math.round(html.length / 1024));
