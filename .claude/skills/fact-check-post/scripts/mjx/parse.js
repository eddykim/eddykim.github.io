// stdin: JSON [{k, tex, display}] → stdout: JSON [{k, error}] (오류만)
const {mathjax} = require('mathjax-full/js/mathjax.js');
const {TeX} = require('mathjax-full/js/input/tex.js');
const {liteAdaptor} = require('mathjax-full/js/adaptors/liteAdaptor.js');
const {RegisterHTMLHandler} = require('mathjax-full/js/handlers/html.js');
const {AllPackages} = require('mathjax-full/js/input/tex/AllPackages.js');
const {SerializedMmlVisitor} = require('mathjax-full/js/core/MmlTree/SerializedMmlVisitor.js');

const {STATE} = require('mathjax-full/js/core/MathItem.js');
const adaptor = liteAdaptor();
RegisterHTMLHandler(adaptor);
let last = null;
const tex = new TeX({
  packages: AllPackages.filter(p => p !== 'noundefined' && p !== 'bussproofs'),
  tags: 'ams',
  formatError: (jax, err) => { last = err.message; throw err; },
});
const doc = mathjax.document('', {InputJax: tex});
const visitor = new SerializedMmlVisitor();
const items = JSON.parse(require('fs').readFileSync(0, 'utf8'));
const out = [];
for (const it of items) {
  last = null;
  try {
    const node = doc.convert(it.tex, {display: it.display, end: STATE.COMPILED});
    const mml = visitor.visitTree(node);
    if (/<merror/.test(mml)) out.push({k: it.k, error: 'merror'});
  } catch (e) {
    out.push({k: it.k, error: last || String(e.message || e)});
  }
}
process.stdout.write(JSON.stringify(out));
