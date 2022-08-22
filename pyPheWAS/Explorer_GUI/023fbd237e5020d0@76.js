// https://observablehq.com/@ddspog/notebook-itens@76
import define1 from "./29d10840f83e2527@465.js";

function _1(md){return(
md`# Notebook Itens`
)}

function _2(cheery){return(
cheery
)}

function _3(md){return(
md`_This notebook export some functions to print, display interesting and useful things onto notebook._`
)}

function _4(md){return(
md`### banner (url, height)
* **url**: an url to an image as string.
* **height**: the height to banner as string.

This function will create an banner with full width on the notebook, resizing image to fit banner, filling its width with repeat.
`
)}

function _banner(html){return(
(url, height, repeat) => {
  if (repeat === undefined) {
    repeat = "repeat-x";
  }
  return html`<div style='height: ${height}; background-image: url(${url}); background-size: auto 256px; background-repeat: ${repeat}; background-position: center;'></div>`;
}
)}

function _6(md){return(
md`### placeholder (width, height, text)
* **width**: the width to placeholder.
* **height**: the height to placeholder.
* **text**: the text to put on placeholder.

This function will create an placeholder image painted on gray, with an optional text, or the dimensions printed at the center.
`
)}

function _placeholder(html,isTextValid){return(
(width, height, text) => {
  return html`<img src="https://via.placeholder.com/${width}x${height}.png${
  isTextValid(text) ? "?text=" + text : ""
  }">`
}
)}

function _8(md){return(
md`---

# Appendix

_Importing Cheery graphics to use as cover._`
)}

function _10(md){return(
md`_Some short functions._`
)}

function _isTextValid(){return(
(text) => {
  return text !== undefined && text !== null && text != ""
}
)}

export default function define(runtime, observer) {
  const main = runtime.module();
  main.variable(observer()).define(["md"], _1);
  main.variable(observer()).define(["cheery"], _2);
  main.variable(observer()).define(["md"], _3);
  main.variable(observer()).define(["md"], _4);
  main.variable(observer("banner")).define("banner", ["html"], _banner);
  main.variable(observer()).define(["md"], _6);
  main.variable(observer("placeholder")).define("placeholder", ["html","isTextValid"], _placeholder);
  main.variable(observer()).define(["md"], _8);
  const child1 = runtime.module(define1);
  main.import("svg", "cheery", child1);
  main.variable(observer()).define(["md"], _10);
  main.variable(observer("isTextValid")).define("isTextValid", _isTextValid);
  return main;
}
