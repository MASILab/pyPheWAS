// https://observablehq.com/@observablehq/cheery-logo@465
export default function define(runtime, observer) {
  const main = runtime.module();
  main.variable(observer()).define(["md"], function(md){return(
md`# Cheery Logo`
)});
  main.variable(observer("svg")).define("svg", ["d3","DOM","width","height","circles","betterPastel","padRadius","icon"], function(d3,DOM,width,height,circles,betterPastel,padRadius,icon)
{
  const svg = d3
    .select(DOM.svg(width, height))
    .attr("viewBox", `${-width / 2} ${-height / 2} ${width} ${height}`)
    .style("width", "100%")
    .style("height", "auto")
    .style("display", "block");

  const circle = svg
    .selectAll("g")
    .data(circles)
    .join("g")
    .attr("transform", d => `translate(${d.x},${d.y})`);

  circle
    .append("circle")
    .attr("fill", (d, i) => (d.logo ? "none" : betterPastel[i % 3]))
    .attr("r", 0)
    .on("click", function() {
      d3.select(this)
        .transition()
        .duration(100)
        .ease(d3.easeCircleIn)
        .attr("r", 0);
    })
    .transition()
    .ease(d3.easeExpIn)
    .delay((d, i) => i * 3)
    .attr("r", d => Math.max(1, d.r - padRadius));

  circle
    .filter(d => d.logo)
    .attr("transform", "translate(-120, -120)")
    .append("path")
    .attr("fill", "#999")
    .attr("transform", d => `scale(${(d.r - padRadius) / 12})`)
    .attr("d", icon);

  return svg.node();
}
);
  main.variable(observer("betterPastel")).define("betterPastel", function(){return(
["#b0deff", "#ffd19a", "#fff8a6"]
)});
  main.variable(observer("padRadius")).define("padRadius", function(){return(
5
)});
  main.variable(observer("height")).define("height", function(){return(
320
)});
  main.variable(observer("scale")).define("scale", function(){return(
15
)});
  main.variable(observer("circles")).define("circles", ["d3","scale"], function(d3,scale)
{
  let circles = d3.packSiblings(
    [{ r: 125, logo: true }].concat(
      d3
        .range(550)
        .map(i => ({ r: 10 + (Math.random() * i) / scale }))
        .reverse()
    )
  );

  const center = circles.find(c => c.logo);

  for (let c of circles) {
    if (c.logo) continue;
    c.x -= center.x;
    c.y -= center.y;
  }

  center.x = center.y = 0;

  return circles;
}
);
  main.variable(observer("icon")).define("icon", function(){return(
`M12.1450213,20.7196596 C11.0175263,20.7196596 10.0411956,20.4623004 9.216,19.9475745 C8.39080438,19.4328485 7.75761923,18.7343023 7.31642553,17.8519149 C6.87523184,16.9695275 6.55251166,16.0340475 6.34825532,15.0454468 C6.14399898,14.0568461 6.04187234,12.990644 6.04187234,11.8468085 C6.04187234,10.9971021 6.09497819,10.1841741 6.20119149,9.408 C6.30740479,8.63182591 6.50348793,7.84340826 6.78944681,7.0427234 C7.07540569,6.24203855 7.44306158,5.54757741 7.89242553,4.95931915 C8.34178948,4.37106089 8.93003892,3.89310822 9.65719149,3.52544681 C10.3843441,3.1577854 11.2136124,2.97395745 12.1450213,2.97395745 C13.2725163,2.97395745 14.2488469,3.23131658 15.0740426,3.74604255 C15.8992382,4.26076853 16.5324233,4.95931474 16.973617,5.84170213 C17.4148107,6.72408952 17.7375309,7.65956953 17.9417872,8.64817021 C18.1460436,9.6367709 18.2481702,10.702973 18.2481702,11.8468085 C18.2481702,12.6965149 18.1950644,13.5094429 18.0888511,14.285617 C17.9826378,15.0617911 17.7824696,15.8502088 17.4883404,16.6508936 C17.1942113,17.4515785 16.8265554,18.1460396 16.3853617,18.7342979 C15.944168,19.3225561 15.3600036,19.8005088 14.6328511,20.1681702 C13.9056985,20.5358316 13.0764302,20.7196596 12.1450213,20.7196596 Z M14.245196,13.9469832 C14.8285807,13.3635984 15.1202688,12.6635472 15.1202688,11.8468085 C15.1202688,11.0300698 14.8358729,10.3300186 14.2670728,9.74663382 C13.6982726,9.16324904 12.9909292,8.87156103 12.1450213,8.87156103 C11.2991134,8.87156103 10.5917699,9.16324904 10.0229698,9.74663382 C9.45416961,10.3300186 9.1697738,11.0300698 9.1697738,11.8468085 C9.1697738,12.6635472 9.45416961,13.3635984 10.0229698,13.9469832 C10.5917699,14.530368 11.2991134,14.822056 12.1450213,14.822056 C12.9909292,14.822056 13.6909804,14.530368 14.245196,13.9469832 Z M12,24 C18.627417,24 24,18.627417 24,12 C24,5.372583 18.627417,0 12,0 C5.372583,0 0,5.372583 0,12 C0,18.627417 5.372583,24 12,24 Z`
)});
  main.variable(observer("d3")).define("d3", ["require"], function(require){return(
require("d3@5")
)});
  return main;
}
