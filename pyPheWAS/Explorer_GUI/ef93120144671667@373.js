// https://observablehq.com/@ddspog/useful-libs@373
import define1 from "./023fbd237e5020d0@76.js";

function _1(md,banner){return(
md`# Useful Libs

${banner("https://cdn-images-1.medium.com/max/1600/1*aalTd6nKuVR31c3-bEED8g.jpeg", "256px")}

_This notebook is exporting libraries used on others Notebooks, in order to avoid using extra blocks to import recorrent libraries._`
)}

function _2(md){return(
md`## Vega-Lite Embed

I've imported this Vega-Lite embed function. The API is described [here](https://github.com/vega/vega-embed). You can import using:

\`\`\`js
import { VegaEmbed as embed } from "@ddspog/useful-libs"
\`\`\`

`
)}

function _VegaEmbed(require){return(
require("vega-embed@3")
)}

function _4(md){return(
md`## JQuery

Simple and old jQuery.

\`\`\`js
import { jQuery as $ } from "@ddspog/useful-libs"
\`\`\`

`
)}

function _jQuery(require){return(
require("jquery")
)}

function _6(md){return(
md`## D3

The famous D3 Library.

\`\`\`js
import { d3 } from "@ddspog/useful-libs"
\`\`\``
)}

async function _d3(require)
{
  var d3 = await require("d3");

  d3.transition.prototype.at = d3.selection.prototype.at;
  d3.transition.prototype.st = d3.selection.prototype.st;
  d3.transition.prototype.tspans = d3.selection.prototype.tspans;
  
  d3.selection.prototype.apply = function(fn, ...attrs) {
    return fn(this, ...attrs);
  }
  
  d3.fetch = (d3Fn, url, row) => new Promise(resolve => {
    let fs = d3Fn(url);
    
    if (row !== undefined)
      fs.row(row);
      
    return fs.get(d => resolve(d))
  });
  
  return d3;
}


function _8(md){return(
md`## Luxon Tools

Utility functions.

\`\`\`js
import { DateTime, Interval } from "@ddspog/useful-libs"
\`\`\``
)}

function _luxon(){return(
import('https://unpkg.com/luxon@2.0.2/src/luxon.js?module')
)}

function _DateTime(luxon){return(
luxon.DateTime
)}

function _Interval(luxon){return(
luxon.Interval
)}

function _12(md){return(
md`## Lodash

Utility functions.

\`\`\`js
import { _ } from "@ddspog/useful-libs"
\`\`\`

`
)}

function __(require){return(
require('lodash')
)}

function _14(md){return(
md`## Math

_Package for simple math operations._

\`\`\`js
import { math } from "@ddspog/useful-libs"
\`\`\`

`
)}

function _math(require){return(
require("https://unpkg.com/mathjs@5.2.0/dist/math.min.js")
)}

function _16(md){return(
md`## Other functions

_Useful functions such as Clone, to safely clone any object, array or similar element._

\`\`\`js
import { Clone } from "@ddspog/useful-libs"
\`\`\``
)}

function _Clone(require){return(
require('https://bundle.run/clone@2.1.2')
)}

function _18(md){return(
md`
_The set function it's a utility to set nested properties values with ease._

\`\`\`js
import { set } from "@ddspog/useful-libs"
\`\`\``
)}

function _set(require){return(
require('https://bundle.run/set-value@3.0.0')
)}

function _20(md){return(
md`

_The get function, similarly to set, access nested properties values._

\`\`\`js
import { get } from "@ddspog/useful-libs"
\`\`\``
)}

function _get(require){return(
require('https://bundle.run/get-value@3.0.1')
)}

function _22(md){return(
md`---
# Appendix

`
)}

function _23(md){return(
md`_Importing banner function to use as cover._`
)}

export default function define(runtime, observer) {
  const main = runtime.module();
  main.variable(observer()).define(["md","banner"], _1);
  main.variable(observer()).define(["md"], _2);
  main.variable(observer("VegaEmbed")).define("VegaEmbed", ["require"], _VegaEmbed);
  main.variable(observer()).define(["md"], _4);
  main.variable(observer("jQuery")).define("jQuery", ["require"], _jQuery);
  main.variable(observer()).define(["md"], _6);
  main.variable(observer("d3")).define("d3", ["require"], _d3);
  main.variable(observer()).define(["md"], _8);
  main.variable(observer("luxon")).define("luxon", _luxon);
  main.variable(observer("DateTime")).define("DateTime", ["luxon"], _DateTime);
  main.variable(observer("Interval")).define("Interval", ["luxon"], _Interval);
  main.variable(observer()).define(["md"], _12);
  main.variable(observer("_")).define("_", ["require"], __);
  main.variable(observer()).define(["md"], _14);
  main.variable(observer("math")).define("math", ["require"], _math);
  main.variable(observer()).define(["md"], _16);
  main.variable(observer("Clone")).define("Clone", ["require"], _Clone);
  main.variable(observer()).define(["md"], _18);
  main.variable(observer("set")).define("set", ["require"], _set);
  main.variable(observer()).define(["md"], _20);
  main.variable(observer("get")).define("get", ["require"], _get);
  main.variable(observer()).define(["md"], _22);
  main.variable(observer()).define(["md"], _23);
  const child1 = runtime.module(define1);
  main.import("banner", child1);
  return main;
}
