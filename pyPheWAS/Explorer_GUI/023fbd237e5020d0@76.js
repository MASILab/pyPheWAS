// https://observablehq.com/@ddspog/notebook-itens@76
import define1 from "./29d10840f83e2527@465.js";

export default function define(runtime, observer) {
  const main = runtime.module();
  main.variable(observer()).define(["md"], function(md){return(
md`# Notebook Itens`
)});
  main.variable(observer()).define(["cheery"], function(cheery){return(
cheery
)});
  main.variable(observer()).define(["md"], function(md){return(
md`_This notebook export some functions to print, display interesting and useful things onto notebook._`
)});
  main.variable(observer()).define(["md"], function(md){return(
md`### banner (url, height)
* **url**: an url to an image as string.
* **height**: the height to banner as string.

This function will create an banner with full width on the notebook, resizing image to fit banner, filling its width with repeat.
`
)});
  main.variable(observer("banner")).define("banner", ["html"], function(html){return(
(url, height, repeat) => {
  if (repeat === undefined) {
    repeat = "repeat-x";
  }
  return html`<div style='height: ${height}; background-image: url(${url}); background-size: auto 256px; background-repeat: ${repeat}; background-position: center;'></div>`;
}
)});
  main.variable(observer()).define(["md"], function(md){return(
md`### placeholder (width, height, text)
* **width**: the width to placeholder.
* **height**: the height to placeholder.
* **text**: the text to put on placeholder.

This function will create an placeholder image painted on gray, with an optional text, or the dimensions printed at the center.
`
)});
  main.variable(observer("placeholder")).define("placeholder", ["html","isTextValid"], function(html,isTextValid){return(
(width, height, text) => {
  return html`<img src="https://via.placeholder.com/${width}x${height}.png${
  isTextValid(text) ? "?text=" + text : ""
  }">`
}
)});
  main.variable(observer()).define(["md"], function(md){return(
md`---

# Appendix

_Importing Cheery graphics to use as cover._`
)});
  const child1 = runtime.module(define1);
  main.import("svg", "cheery", child1);
  main.variable(observer()).define(["md"], function(md){return(
md`_Some short functions._`
)});
  main.variable(observer("isTextValid")).define("isTextValid", function(){return(
(text) => {
  return text !== undefined && text !== null && text != ""
}
)});
  return main;
}
