// https://observablehq.com/@ddspog/useful-libs@352
import define1 from "./023fbd237e5020d0@76.js";

export default function define(runtime, observer) {
  const main = runtime.module();
  main.variable(observer()).define(["md","banner"], function(md,banner){return(
md`# Useful Libs

${banner("https://cdn-images-1.medium.com/max/1600/1*aalTd6nKuVR31c3-bEED8g.jpeg", "256px")}

_This notebook is exporting libraries used on others Notebooks, in order to avoid using extra blocks to import recorrent libraries._`
)});
  main.variable(observer()).define(["md"], function(md){return(
md`## Vega-Lite Embed

I've imported this Vega-Lite embed function. The API is described [here](https://github.com/vega/vega-embed). You can import using:

\`\`\`js
import { VegaEmbed as embed } from "@ddspog/useful-libs"
\`\`\`

`
)});
  main.variable(observer("VegaEmbed")).define("VegaEmbed", ["require"], function(require){return(
require("vega-embed@3")
)});
  main.variable(observer()).define(["md"], function(md){return(
md`## JQuery

Simple and old jQuery.

\`\`\`js
import { jQuery as $ } from "@ddspog/useful-libs"
\`\`\`

`
)});
  main.variable(observer("jQuery")).define("jQuery", ["require"], function(require){return(
require("jquery")
)});
  main.variable(observer()).define(["md"], function(md){return(
md`## D3

The famous D3 Library.

\`\`\`js
import { d3 } from "@ddspog/useful-libs"
\`\`\``
)});
  main.variable(observer("d3")).define("d3", ["require"], async function(require)
{
  var d3 = await require("d3-jetpack/build/d3v4+jetpack");

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
);
  main.variable(observer()).define(["md"], function(md){return(
md`## Lodash

Utility functions.

\`\`\`js
import { _ } from "@ddspog/useful-libs"
\`\`\`

`
)});
  main.variable(observer("_")).define("_", ["require"], function(require){return(
require('lodash')
)});
  main.variable(observer()).define(["md"], function(md){return(
md`## Math

_Package for simple math operations._

\`\`\`js
import { math } from "@ddspog/useful-libs"
\`\`\`

`
)});
  main.variable(observer("math")).define("math", ["require"], function(require){return(
require("https://unpkg.com/mathjs@5.2.0/dist/math.min.js")
)});
  main.variable(observer()).define(["md"], function(md){return(
md`## Other functions

_Useful functions such as Clone, to safely clone any object, array or similar element._

\`\`\`js
import { Clone } from "@ddspog/useful-libs"
\`\`\``
)});
  main.variable(observer("Clone")).define("Clone", ["require"], function(require){return(
require('https://bundle.run/clone@2.1.2')
)});
  main.variable(observer()).define(["md"], function(md){return(
md`
_The set function it's a utility to set nested properties values with ease._

\`\`\`js
import { set } from "@ddspog/useful-libs"
\`\`\``
)});
  main.variable(observer("set")).define("set", ["require"], function(require){return(
require('https://bundle.run/set-value@3.0.0')
)});
  main.variable(observer()).define(["md"], function(md){return(
md`

_The get function, similarly to set, access nested properties values._

\`\`\`js
import { get } from "@ddspog/useful-libs"
\`\`\``
)});
  main.variable(observer("get")).define("get", ["require"], function(require){return(
require('https://bundle.run/get-value@3.0.1')
)});
  main.variable(observer()).define(["md"], function(md){return(
md`---
# Appendix

`
)});
  main.variable(observer()).define(["md"], function(md){return(
md`_Importing banner function to use as cover._`
)});
  const child1 = runtime.module(define1);
  main.import("banner", child1);
  return main;
}
