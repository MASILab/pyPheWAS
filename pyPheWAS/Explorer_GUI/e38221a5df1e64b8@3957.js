import define1 from "./ef93120144671667@373.js";

function _1(md){return(
md `# pyPheWAS Explorer
#### A visualization tool for exploratory analysis of phenome-disease association 
`
)}

function _2(group_panel_size,margin,d3,group_panel)
{
  // Master SVG Element for the Regression Builder Panel
  let full_width = group_panel_size+(2*margin.y)
  let full_height = group_panel_size+(2*margin.y)
  let svg_elem = d3.create('svg')
    .attr('id', 'group_panel')
    .attr('width', full_width)
    .attr('height', full_height)
    .style('font-family', 'sans-serif')
    .style('font-size', 10)
  
  let reg_group = svg_elem.append('g').attr('transform', 'translate('+(group_panel_size+margin.x)+','+margin.y+')')
  
  group_panel (svg_elem)
  
  return svg_elem.node()
  }


function _3(reg_panel_size,margin,d3,reg_panel)
{
  // Master SVG Element for the Regression Evaluation Panel
  let full_width = reg_panel_size+(2*margin.y)
  let full_height = reg_panel_size+(2*margin.y)
  let svg_elem = d3.create('svg')
    .attr('id', 'reg_panel')
    .attr('width', full_width)
    .attr('height', full_height)
    .style('font-family', 'sans-serif')
    .style('font-size', 10)

  reg_panel(svg_elem)
  
  return svg_elem.node()
  }


function _4(md){return(
md `### Group Panel Functions`
)}

function _group_panel(group_panel_size,margin,var_dist_size,legend_height,plot_groupvar_hists,indep_comp,reg_builder,group_legend,d3,var_dist2_size){return(
function group_panel (svg_elem) {
  // Setup Plot
  let full_width = group_panel_size+(2*margin.x)
  let full_height = group_panel_size+(2*margin.y)
  
  let main_group = svg_elem.append('g')
    .attr('transform', 'translate('+margin.x+','+(margin.y)+')')
    .attr('id', 'main')

  let var_dist_group = main_group.append('g')
    .attr('id', 'var_dist_panel') 

  let legend_group = main_group.append('g')
    .attr('transform', 'translate('+(var_dist_size.width)+',0)')
    .attr('id', 'legend') 

  let var_comp_group = main_group.append('g')
    .attr('transform', 'translate('+(var_dist_size.width)+','+(legend_height)+')')
    .attr('id', 'var_comp_panel')
 
  plot_groupvar_hists(var_dist_group)
  indep_comp(var_comp_group)
  reg_builder(main_group)
  group_legend(legend_group)
  
  var_dist_group.append('rect')
    .attr('width', var_dist_size.width)
    .attr('height', var_dist_size.height)
    .attr('stroke', d3.hcl(0,0,70))
    .attr('stroke-width', 2)
    .attr('fill','none')
    .lower()

  var_comp_group.append('rect')
    .attr('width', var_dist2_size.width)
    .attr('height', var_dist2_size.height)
    .attr('stroke', d3.hcl(0,0,70))
    .attr('stroke-width', 2)
    .attr('fill','none')
    .lower()
    
  return svg_elem.node()
}
)}

function _plot_groupvar_hists(d3,group_vars,var_dist_size,z,hist_data,group_data,target_var,$0,$1,$2,cov_select,$3,var_comp,$4,var_stat_scale,draw_groupvar_buttons){return(
function plot_groupvar_hists(hist_group) {
  
  let hgroups = []
  let hgroup_scale = d3.scaleBand().domain(group_vars).range([0,var_dist_size.height])
  
  const gvar_buffer = 5;
  const cell_height = hgroup_scale(group_vars[1]) - hgroup_scale(group_vars[0]);
  const cell_width = var_dist_size.width;
  const hist_height = cell_height / 2.3;
  const hist_width = 100;
  const hist_margin = 37;
  const cor_height = cell_height / 3.0;
  const cor_width = 30;
  var y_offset1 = (cell_height/2) - (hist_height / 2) - (gvar_buffer*2); // approximately 1/4 of the way down the cell
  var y_offset2 = (cell_height/2) - (gvar_buffer*2); // approximately 1/2 of the way down the cell
  var box_offset_x = 10 + (cell_width * (1/2)) // approximately 1/2 of the way across the cell
  let button_size = 11;
  
  for (var i=0; i < group_vars.length; i++){
    let v = group_vars[i]
    let v_data = z.filter(d => d.var_name === v, hist_data)
    let v_data0 = z.filter(d => d.response === 0, v_data)
    let v_data1 = z.filter(d => d.response === 1, v_data)
    
    hgroups[i] = hist_group.append('g')
      .attr("transform", `translate(${gvar_buffer},${hgroup_scale(v)+gvar_buffer})`)
      .attr('id', v)
    
    // seperating line
    if (i !== 0) {
      hgroups[i].append('line')
        .attr('x1',-gvar_buffer)
        .attr('x2',cell_width-gvar_buffer)
        .attr('y1',-gvar_buffer)
        .attr('y2',-gvar_buffer)
        .attr('stroke', d3.hcl(0,0,70))
        .attr('stroke-width',2)
    }
    
    hgroups[i].append('text')
      .text(v)
      .attr('class', v)
      .attr('x', 10+(cell_width*(1/2)))
      .attr('y', y_offset1 + gvar_buffer)
      .attr('font-size', '12px')
      .attr('font-weight', 'bold')
      .attr('text-anchor','start')

    if (group_data[i].ropt) {
      let x_delta = cell_width*0.94
      let button_color = 'gray';
      if (v == target_var){button_color = 'yellow'}
      // Target Variable Button
      hgroups[i].append('rect')
        .attr('class','target_button')
        .attr('transform', `translate(${x_delta - gvar_buffer},${y_offset2})`)
        .attr('width',button_size)
        .attr('height',button_size)
        .attr('stroke',d3.hcl(0,0,70))
        .attr('stroke-width', 1)
        .attr("fill", button_color)
        .attr("opacity",0.7)
        .on('click', function(d) {
          if (v !== target_var){
            $0.value = v
            $1.value = [{cmd:"group_data", response: v}];
            $2.value = [{cmd:"histograms", response: v}];
            // if currently selected as covariate - deselect
            // preserve if any other variables are selected
            let curr_cov_select = z.getCol("s", cov_select)
            let vc_id = z.getCol("id", z.filter(r => r.vname === v, cov_select))[0]
            if (curr_cov_select[vc_id] === 1){
              curr_cov_select[vc_id] = 0
              $3.value = z.addCol("s",curr_cov_select,cov_select)
            }
            // if currently selected for comparison - deselect
            // preserve if any other variables are selected
            let curr_select = z.getCol("s",var_comp)
            let v_id = z.getCol("id", z.filter(r => r.vname === v, var_comp))[0]
            if (curr_select[v_id] === 1){
              let num_selected = z.filter(r => r.s===1,var_comp).length
              if (num_selected === 2){
                let s_id = z.filter(r => r.s===1,var_comp)[0].id
                if (s_id === v_id) { s_id = z.filter(r => r.s===1,var_comp)[1].id }
                curr_select = new Array(curr_select.length).fill(0)
                curr_select[s_id] = 1
              }
              else {curr_select[v_id] = 0}
              $4.value = z.addCol("s",curr_select,var_comp)
            }
          }
      })
      hgroups[i].append('text')
          .text("Target")
          .attr('transform', d=> `translate(${x_delta - gvar_buffer - 3},${y_offset2 + gvar_buffer*2})`)
          .attr('font-size', '11px')
          .attr('text-anchor', 'end')
    }

    // target variable correlation
    hgroups[i].append('rect')
      .attr('class', v)
      .attr('x', box_offset_x)
      .attr('y', y_offset2)
      .attr('width', cor_width)
      .attr('height', cor_height)
      .attr('fill', var_stat_scale(group_data[i].pval))
      .attr('stroke', d3.hcl(0,0,70))
      .on('mouseover', function(d) {
        let this_group = d3.select('#'+v);
        this_group.append('text')
          .attr("class", "corr_description")
          .text("Correlation with "+target_var) 
          .attr('text-anchor', 'start')
          .attr('x', +d3.select(this).attr('x'))
          .attr('y', +d3.select(this).attr('y')-5)
          .attr('font-size', '11px')
          .attr('fill', 'darkslategrey')
          .raise()
      })
      .on('mouseout', (d, i, arr) => {
        let this_group = d3.select('#'+v);
        this_group.selectAll('.corr_description').remove()
      })

    hgroups[i].append('text')
      .text(group_data[i].corr.toFixed(2))
      .attr('font-size','12px')
      .attr('font-weight','bold')
      .attr('text-anchor','middle')
      .attr('x', box_offset_x + cor_width*0.5)
      .attr('y', y_offset2 + cor_height*0.6)

    var x_scale = d3.scaleLinear()
        .domain([d3.min(v_data, d => d.xmin),d3.max(v_data, d => d.xmax)]).nice()
        .range([0, hist_width]);

    var y_scale = d3.scaleLinear()
      .domain([0, d3.max(v_data, d => d.count)]).nice()
      .range([hist_height,0])
    
    hgroups[i].append("g")
      .attr('class', v)
      .attr("transform", `translate(${hist_margin},${y_offset1})`)
      .attr("fill", "purple")
      .attr("opacity",0.6)
      .selectAll("empty").data(v_data0).enter().append("rect")
        .attr("x", d => x_scale(d.xmin))
        .attr("y", d=> y_scale(d.count))
        .attr("width", d => x_scale(d.xmax) - x_scale(d.xmin))
        .attr("height", d=> hist_height - y_scale(d.count))
    
    hgroups[i].append("g")
      .attr('class', v)
      .attr("transform", `translate(${hist_margin},${y_offset1})`)
      .attr("fill", "green")
      .attr("opacity",0.4)
      .selectAll("empty").data(v_data1).enter().append("rect")
        .attr("x", d => x_scale(d.xmin))
        .attr("y", d=> y_scale(d.count))
        .attr("width", d => x_scale(d.xmax) - x_scale(d.xmin))
        .attr("height",d=> hist_height - y_scale(d.count))
    
    const yAxis = g => g
      .attr("transform", `translate(${hist_margin},${y_offset1})`)
      .call(d3.axisLeft(y_scale).ticks(3))
      .call(g => g.select(".domain").remove())
    const xAxis = g => g
      .attr("transform", `translate(${hist_margin},${y_offset1 + hist_height})`)
      .call(d3.axisBottom(x_scale).ticks(5).tickSizeOuter(0))

    hgroups[i].append("g")
      .attr('class', v)
      .call(xAxis);

    hgroups[i].append("g")
      .attr('class', v)
      .call(yAxis);
  }
  
  draw_groupvar_buttons(hist_group, hgroup_scale, cell_height, cell_width, y_offset2+(2*gvar_buffer)+button_size, button_size)
}
)}

function _draw_groupvar_buttons(var_comp,d3,target_var,z,$0,cov_select,$1){return(
function draw_groupvar_buttons(hist_group, groupvar_scale, cell_height, cell_width, y_offset, button_size) {
  let button_color = ["grey","yellow","black"]
  
 // Variable Comparison Buttons
  let varcomp_button_group = hist_group.append('g')
    .selectAll("empty").data(var_comp).enter()

  varcomp_button_group.append('rect')
    .attr('class','varcomp_button')
    .attr('transform', d=> `translate(${cell_width*0.8},${groupvar_scale(d.vname) + y_offset})`)
    .attr('width',button_size)
    .attr('height',button_size)
    .attr('stroke',d3.hcl(0,0,70))
    .attr('stroke-width', 1)
    .attr("fill", function(d){
      if (d.vname == target_var){ return button_color[2] }
      return button_color[d.s];
    })
    .attr("opacity",0.7)
    .on('click', function(d) {
      if (d.vname !== target_var){
        let curr_select = z.getCol("s",var_comp)
        let num_selected = z.filter(r => r.s===1,var_comp).length
        if (curr_select[d.id] === 1){
          curr_select[d.id] = 0
          if (num_selected === 2){
            let s_id = z.filter(r => r.s===1,var_comp)[0].id
            if (s_id === d.id) { s_id = z.filter(r => r.s===1,var_comp)[1].id }
            curr_select = new Array(curr_select.length).fill(0)
            curr_select[s_id] = 1
          }
          $0.value = z.addCol("s",curr_select,var_comp)
        }
        else if (num_selected === 0) {
          curr_select[d.id] = 1
          $0.value = z.addCol("s",curr_select,var_comp)
        }
        else if (num_selected === 1) {
          let s_id = z.filter(r => r.s===1,var_comp)[0].id
          curr_select = new Array(curr_select.length).fill(2)
          curr_select[d.id] = 1
          curr_select[s_id] = 1
          $0.value = z.addCol("s",curr_select,var_comp)
        }
      }
  })
  
  let x_delta = cell_width*0.8  - 4
    varcomp_button_group.append('text')
      .text("Comp")
      .attr('transform', d=> `translate(${x_delta},${y_offset+10+groupvar_scale(d.vname)})`)
      .attr('font-size', '11px')
      .attr('text-anchor', 'end')
  
  // Covariate Buttons
  let cov_button_group = hist_group.append('g')
    .selectAll("empty").data(cov_select).enter()

  cov_button_group.append('rect')
    .attr('class','varcomp_button')
    .attr('transform', d=> `translate(${cell_width*0.94},${groupvar_scale(d.vname)+y_offset})`)
    .attr('width',button_size)
    .attr('height',button_size)
    .attr('stroke',d3.hcl(0,0,70))
    .attr('stroke-width', 1)
    .attr("fill", function(d){
      if (d.vname == target_var){ return button_color[2] }
      return button_color[d.s];
    })
    .attr("opacity",0.7)
    .on('click', function(d) {
      if (d.vname !== target_var){
        let curr_select = z.getCol("s",cov_select)
        if (curr_select[d.id] === 1){
          curr_select[d.id] = 0
          $1.value = z.addCol("s",curr_select,cov_select)
        }
        else {
          curr_select[d.id] = 1
          $1.value = z.addCol("s",curr_select,cov_select)
        }
      }
  })

  x_delta = cell_width*0.94  - 3
    varcomp_button_group.append('text')
      .text("Cov")
      .attr('transform', d=> `translate(${x_delta},${y_offset+10+groupvar_scale(d.vname)})`)
      .attr('font-size', '11px')
      .attr('text-anchor', 'end')

  
}
)}

function _indep_comp(jhist_data,var_dist2_size,z,var1_select,hist_data,var2_select,d3,comp_stats){return(
function indep_comp(var_comp_group) {
  if (jhist_data.length < 200 && jhist_data[0].msg === "no_data"){
    var_comp_group.append('text')
      .attr('class','placeholder')
      .text('Group Variable Comparison')
      .attr('text-anchor','middle')
      .attr('transform',`translate(${var_dist2_size.width/2},150)`)
      .attr('font-size','24px')
      .attr('font-weight','bold')
  }
  else if (jhist_data.length < 200 && jhist_data[0].msg === "select_2nd_var"){
    var_comp_group.append('text')
      .attr('class','placeholder')
      .text('Please select a second variable for comparison')
      .attr('text-anchor','middle')
      .attr('transform',`translate(${var_dist2_size.width/2},150)`)
      .attr('font-size','18px')
      .attr('font-weight','bold')
  }
  else{
    var_comp_group.selectAll('.placeholder').remove()
    
    let var1_bins = z.unique(z.getCol("v1_min",jhist_data)).length
    let var2_bins = z.unique(z.getCol("v2_min",jhist_data)).length
    let var1_data = z.filter(d => d.var_name === var1_select, hist_data)
    let var2_data = z.filter(d => d.var_name === var2_select, hist_data)

    const hist_height = 25;
    const jhist_size = 120;
    const jhist_cellsize = jhist_size / var1_bins;
    const mw = 150;
    const mh = 10;
    let panel_scale = d3.scaleBand()
      .domain(['jh1','jh2','stats'])
      .range([10,var_dist2_size.height-10])
      .paddingInner(0.1)
    let square_scale = d3.scaleBand()
      .domain(d3.range(var1_bins))
      .range([0,jhist_size])

    // Axis Labels
    var_comp_group.append("text")
      .text(var1_select)
      .attr('text-anchor','middle')
      .attr("transform",`translate(${120},${panel_scale('jh2')}),rotate(270)`)
      .attr('font-size', '14px')
      .attr('font-weight', 'bold')

    var_comp_group.append("text")
      .text(var2_select)
      .attr('text-anchor','middle')
      .attr("transform",`translate(${mw+(jhist_size/2)},${panel_scale('jh2')+jhist_size+hist_height+mh+20})`)
      .attr('font-size', '14px')
      .attr('font-weight', 'bold')

    // Genotype = 0
    let g0_group = var_comp_group.append("g")
      .attr("transform", `translate(${mw},${panel_scale('jh1')})`)

    // Joint plot
    let jhist_data0 = z.filter(d => d.response === 0, jhist_data)
    // scales
    let count0_color_scale = d3.scaleLinear()
      .domain([0,1,d3.max(jhist_data0,d=>d.count)])
      .range(['white',d3.hcl(0,0,90), 'purple'])

    g0_group.selectAll("empty").data(jhist_data0).enter().append("rect")
      .attr('transform', (d,i) => {
        var y_off = Math.floor(i/var1_bins)
        var square_offset_y = y_off*jhist_cellsize
        var x_off = i % var1_bins
        var square_offset_x = square_scale(x_off)
        return `translate(${square_offset_x},${square_offset_y})`
      })
      .attr('height', jhist_cellsize)
      .attr('width', jhist_cellsize)
      .attr("fill", d => count0_color_scale(d.count))
      .attr('stroke', d3.hcl(0,0,70))
      .attr('x',0)
      .attr('y',hist_height)
      .on('mouseover', function(d) {
        g0_group.append("text")
          .attr("class","jhist0")
          .text(d.count)
          .attr('text-anchor', 'start')
          .attr('x', jhist_size+5)
          .attr('y', mh+5)
          .attr('font-weight', 'bold')
          .style('font-size', 13)
          .raise()
      })
      .on('mouseout', (d, i, arr) => {
        g0_group.selectAll('.jhist0').remove()
      })

    // Var 1 Histogram (Y-axis)
    let v1_data0 = z.filter(d => d.response === 0, var1_data)
    var x_scale = d3.scaleLinear()
      .domain([d3.min(var1_data, d => d.xmin),d3.max(var1_data, d => d.xmax)])
      .range([0, jhist_size]);
    var y_scale = d3.scaleLinear()
      .domain([0, d3.max(var1_data, d => d.count)]).nice()
      .range([hist_height,0])
    g0_group.append("g")
      .attr("fill", "purple")
      .attr("opacity",0.6)
      .attr("transform",`translate(${jhist_size+hist_height},${hist_height}), rotate(90)`)
      .selectAll("empty").data(v1_data0).enter().append("rect")
        .attr("x", d => x_scale(d.xmin))
        .attr("y", d=> y_scale(d.count))
        .attr("width", d => x_scale(d.xmax) - x_scale(d.xmin))
        .attr("height",d=> hist_height - y_scale(d.count))
    //axis
    const var1_Axis0 = g => g
      .attr("transform", `translate(0,${hist_height}), rotate(90)`)
      .call(d3.axisBottom(x_scale).ticks(5).tickSizeOuter(0))
      .selectAll("text")
        .attr("transform", "rotate(-90)")
        .attr("x",-15)
        .attr("y",-5)
        .style("text-anchor", "middle");
    g0_group.append("g")
      .call(var1_Axis0);

    // Var 2 Histogram (X-axis)
    let v2_data0 = z.filter(d => d.response === 0, var2_data)
    var x_scale = d3.scaleLinear()
      .domain([d3.min(var2_data, d => d.xmin),d3.max(var2_data, d => d.xmax)])
      .range([0, jhist_size]);
    var y_scale = d3.scaleLinear()
      .domain([0, d3.max(var2_data, d => d.count)]).nice()
      .range([hist_height,0])
    g0_group.append("g")
      .attr("fill", "purple")
      .attr("opacity",0.6)
      .selectAll("empty").data(v2_data0).enter().append("rect")
        .attr("x", d => x_scale(d.xmin))
        .attr("y", d=> y_scale(d.count))
        .attr("width", d => x_scale(d.xmax) - x_scale(d.xmin))
        .attr("height",d=> hist_height - y_scale(d.count))


    // Genotype = 1 
    let g1_group = var_comp_group.append("g")
      .attr("transform", `translate(${mw},${panel_scale('jh2')})`)

    // Joint Plot
    let jhist_data1 = z.filter(d => d.response === 1, jhist_data)
    // scales
    let count1_color_scale = d3.scaleLinear()
      .domain([0,1,d3.max(jhist_data1,d=>d.count)])
      .range(['white',d3.hcl(0,0,90), 'green'])

    g1_group.selectAll("empty").data(jhist_data1).enter().append("rect")
      .attr('transform', (d,i) => {
        var y_off = Math.floor(i/var2_bins)
        var square_offset_y = y_off*jhist_cellsize
        var x_off = i % var2_bins
        var square_offset_x = square_scale(x_off)
        return 'translate('+square_offset_x+','+square_offset_y+')'
      })
      .attr('x',0)
      .attr('y',hist_height)
      .attr('height', jhist_cellsize)
      .attr('width', jhist_cellsize)
      .attr("fill", d => count1_color_scale(d.count))
      .attr('stroke', d3.hcl(0,0,70))
      .on('mouseover', function(d) {
        g1_group.append("text")
          .attr("class","jhist1")
          .text(d.count)
          .attr('text-anchor', 'start')
          .attr('x', jhist_size+5)
          .attr('y', mh+5)
          .attr('font-weight', 'bold')
          .style('font-size', 13)
          .raise()
      })
      .on('mouseout', (d, i, arr) => {
        g1_group.selectAll('.jhist1').remove()
      })

    // Var 1 Histogram (Y-axis)
    let v1_data1 = z.filter(d => d.response === 1, var1_data)
    var x_scale = d3.scaleLinear()
      .domain([d3.min(var1_data, d => d.xmin),d3.max(var1_data, d => d.xmax)])
      .range([0, jhist_size]);
    var y_scale = d3.scaleLinear()
      .domain([0, d3.max(var1_data, d => d.count)]).nice()
      .range([hist_height,0])
    g1_group.append("g")
      .attr("fill", "green")
      .attr("opacity",0.6)
      .attr("transform",`translate(${jhist_size+hist_height},${hist_height}), rotate(90)`)
      .selectAll("empty").data(v1_data1).enter().append("rect")
        .attr("x", d => x_scale(d.xmin))
        .attr("y", d=> y_scale(d.count))
        .attr("width", d => x_scale(d.xmax) - x_scale(d.xmin))
        .attr("height",d=> hist_height - y_scale(d.count))
    //axis
    const var1_Axis1 = g => g
      .attr("transform", `translate(0,${hist_height}), rotate(90)`)
      .call(d3.axisBottom(x_scale).ticks(5).tickSizeOuter(0))
      .selectAll("text")
        .attr("transform", "rotate(-90)")
        .attr("x",-15)
        .attr("y",-5)
        .style("text-anchor", "middle");
    g1_group.append("g")
      .call(var1_Axis1);

    // Var 2 Histogram (X-axis)
    let v2_data1 = z.filter(d => d.response === 1, var2_data)
    var x_scale = d3.scaleLinear()
      .domain([d3.min(var2_data, d => d.xmin),d3.max(var2_data, d => d.xmax)])
      .range([0, jhist_size]);
    var y_scale = d3.scaleLinear()
      .domain([0, d3.max(var2_data, d => d.count)]).nice()
      .range([hist_height,0])
    g1_group.append("g")
      .attr("fill", "green")
      .attr("opacity",0.6)
      .selectAll("empty").data(v2_data1).enter().append("rect")
        .attr("x", d => x_scale(d.xmin))
        .attr("y", d=> y_scale(d.count))
        .attr("width", d => x_scale(d.xmax) - x_scale(d.xmin))
        .attr("height",d=> hist_height - y_scale(d.count))
    //axis
    const var2_Axis1 = g => g
      .attr("transform", `translate(0,${jhist_size+hist_height})`)
      .call(d3.axisBottom(x_scale).ticks(5).tickSizeOuter(0))
    g1_group.append("g")
      .call(var2_Axis1);

    let stats_group = var_comp_group.append("g")
      .attr("transform", `translate(0,${panel_scale('stats')})`)

    comp_stats(stats_group)
  }

}
)}

function _comp_stats(d3,var_comp_stats,var_dist2_size,var1_select,var2_select,z,var_stat_scale,target_var){return(
function comp_stats (stats_comp_group){
  const mh = 40
  const mw = 140
  const text_offset = 130
  
  // scales
  let col_scale = d3.scaleBand()
    .domain(d3.set(var_comp_stats, d=> d.test_name).values())
    .range([mw,var_dist2_size.width-20])
    .paddingInner(0.4)
  
  let row_scale = d3.scaleBand()
    .domain([var1_select,var2_select])
    .range([0,60])
    .paddingInner(0.3)
  
  // variable labels
  stats_comp_group.append("text")
    .text(var1_select)
    .attr('text-anchor','end')
    .attr("transform",`translate(${text_offset},${5+mh+row_scale(var1_select)+row_scale.bandwidth()/2})`)
    .attr('font-size', '12px')
    .attr('font-weight', 'bold')
  
  stats_comp_group.append("text")
    .text(var2_select)
    .attr('text-anchor','end')
    .attr("transform",`translate(${text_offset},${5+mh+row_scale(var2_select)+row_scale.bandwidth()/2})`)
    .attr('font-size', '12px')
    .attr('font-weight', 'bold')
  
  // correlation
  let corr_data = z.filter(r => r.test_name === 'correlation', var_comp_stats)
  
  let corr_group = stats_comp_group.append('g')
    .attr('transform',`translate(${col_scale('correlation')},${mh})`)
    .selectAll("empty").data(corr_data).enter()

  corr_group.append('rect')
    .attr('y',row_scale(var1_select))
    .attr('width', col_scale.bandwidth())
    .attr('height', row_scale.range()[1])
    .attr('stroke-width',2)
    .attr('fill', d => var_stat_scale(d.pval))
    .attr('stroke', d3.hcl(0,0,70))
    .on('mouseover', function(d) {
        stats_comp_group.append('text')
          .attr("class", "corr2_description")
          .text("Correlation between "+var1_select+" and "+var2_select) 
          .attr('text-anchor', 'start')
          .attr('x', +d3.select(this).attr('x')+10)
          .attr('y', mh+row_scale.range()[1]+40)
          .attr('font-size', '11px')
          .attr('fill', 'darkslategrey')
          .raise()
      })
      .on('mouseout', (d, i, arr) => {
        stats_comp_group.selectAll('.corr2_description').remove()
      })

  corr_group.append('text')
    .text(d => d.result.toFixed(4))
    .attr('font-size','12px')
    .attr('font-weight','bold')
    .attr('text-anchor','middle')
    .attr('x', col_scale.bandwidth()/2)
    .attr('y',row_scale.range()[1]-25)
  
  corr_group.append('text')
    .text('Correlation')
    .attr('font-size','12px')
    .attr('text-anchor','middle')
    .attr('x', col_scale.bandwidth()/2)
    .attr('y',row_scale.range()[1]+20)
  
  // univariate regression
  let uni_data = z.filter(r => r.test_name === 'univariate regression', var_comp_stats)
  
  let uni_group = stats_comp_group.append('g')
    .attr('transform',`translate(${col_scale('univariate regression')},${mh})`)
  
  uni_group.selectAll("empty").data(uni_data).enter().append('rect')
    .attr('y',d => row_scale(d.var))
    .attr('width', col_scale.bandwidth())
    .attr('height', row_scale.bandwidth())
    .attr('stroke-width',2)
    .attr('fill', d => var_stat_scale(d.pval))
    .attr('stroke', d3.hcl(0,0,70))
    .on('mouseover', function(d) {
        stats_comp_group.append('text')
          .attr("class", "univariate_description")
          .text("Reg coeff for "+d.var+" [ "+target_var+" ~ "+d.var+" ]") 
          .attr('text-anchor', 'start')
          .attr('x', +d3.select(this).attr('x')+10)
          .attr('y', mh+row_scale.range()[1]+40)
          .attr('font-size', '11px')
          .attr('fill', 'darkslategrey')
          .raise()
      })
    .on('mouseout', (d, i, arr) => {
        stats_comp_group.selectAll('.univariate_description').remove()
      })

  uni_group.selectAll("empty").data(uni_data).enter().append('text')
    .text(d => d.result.toFixed(4))
    .attr('font-size','12px')
    .attr('font-weight','bold')
    .attr('text-anchor','middle')
    .attr('x', col_scale.bandwidth()/2)
    .attr('y',d => row_scale(d.var)+(mh/2))
  
  uni_group.append('text')
    .text('Univariate Reg')
    .attr('font-size','12px')
    .attr('text-anchor','middle')
    .attr('x', col_scale.bandwidth()/2)
    .attr('y',row_scale.range()[1]+20)
  
  // multivariate regression
  let mul_data = z.filter(r => r.test_name === 'multivariate regression', var_comp_stats)
  
  let mul_group = stats_comp_group.append('g')
    .attr('transform',`translate(${col_scale('multivariate regression')},${mh})`)
  
  mul_group.selectAll("empty").data(mul_data).enter().append('rect')
    .attr('y',d => row_scale(d.var))
    .attr('width', col_scale.bandwidth())
    .attr('height', row_scale.bandwidth())
    .attr('stroke-width',2)
    .attr('fill', d => var_stat_scale(d.pval))
    .attr('stroke', d3.hcl(0,0,70))
    .on('mouseover', function(d) {
        stats_comp_group.append('text')
          .attr("class", "multivariate_description")
          .text("Reg coeff for "+d.var+" [ "+target_var+" ~ "+var1_select+" + "+var2_select+" ]") 
          .attr('text-anchor', 'start')
          .attr('x', +d3.select(this).attr('x')+10)
          .attr('y', mh+row_scale.range()[1]+40)
          .attr('font-size', '11px')
          .attr('fill', 'darkslategrey')
          .raise()
      })
    .on('mouseout', (d, i, arr) => {
        stats_comp_group.selectAll('.multivariate_description').remove()
      })

  mul_group.selectAll("empty").data(mul_data).enter().append('text')
    .text(d => d.result.toFixed(4))
    .attr('font-size','12px')
    .attr('font-weight','bold')
    .attr('text-anchor','middle')
    .attr('x', col_scale.bandwidth()/2)
    .attr('y',d => row_scale(d.var)+(mh/2))
  
  mul_group.append('text')
    .text('Multivariate Reg')
    .attr('font-size','12px')
    .attr('text-anchor','middle')
    .attr('x', col_scale.bandwidth()/2)
    .attr('y',row_scale.range()[1]+20)
  
}
)}

function _reg_builder(group_data,reg_builder_size,var_dist_size,d3,math,reg_types,$0,z,cov_selection){return(
function reg_builder(main_group){
  
  const num_g0 = group_data[0].g0
  const num_g1 = group_data[0].g1
  const cohort_size = num_g0+num_g1
  
  const geno_offset = {x:30,y:reg_builder_size.height/3 - 5}
  const hist_width = 75
  const geno_panel_width = reg_builder_size.width/4
  const reg_panel_width = reg_builder_size.width*(3/4)
  
  let reg_builder_group = main_group.append('g')
    .attr('transform', 'translate(0,'+var_dist_size.height+')')
    .attr('id', 'reg_builder')
  
 
  // Cohort info
  reg_builder_group.append('rect')
    .attr('width', geno_panel_width)
    .attr('height', reg_builder_size.height)
    .attr('stroke', d3.hcl(0,0,70))
    .attr('stroke-width', 1.5)
    .attr('fill','none')
    .lower()
  reg_builder_group.append('text')
    .text(cohort_size+" patients in cohort")
    .attr('text-anchor', 'middle')
    .attr('x', geno_panel_width/2)
    .attr('y', geno_offset.y -35)
    .attr('font-size', '14px')
    .attr('font-weight', 'bold')
  // reg_builder_group.append('text')
  //   .text(target_var)
  //   .attr('text-anchor', 'middle')
  //   .attr('transform', `rotate(270), translate(${-geno_offset.y},${geno_offset.x})`)
  //   .attr('font-size', '14px')
  let geno_x_scale = d3.scaleLinear()
    .domain([0,math.max(num_g0,num_g1)])
    .range([geno_offset.x+15,geno_offset.x+hist_width+15])
  
  let row_scale = d3.scaleBand()
    .domain([0,1])
    .range([0,50])
    .paddingInner(0.4)
  
  reg_builder_group.append('rect')
    .attr('transform', `translate(0,${geno_offset.y-27+row_scale(0)})`)
    .attr('x',geno_x_scale(0))
    .attr('width',geno_x_scale(num_g0)-geno_x_scale(0))
    .attr('height',row_scale.bandwidth())
    .attr("fill", "purple")
    .attr("opacity",0.6)
  reg_builder_group.append('text')
    .text(num_g0)
    .attr('y',12)
    .attr('transform', `translate(${geno_x_scale(num_g0)+5},${geno_offset.y-27+row_scale(0)})`)
    .attr('font-size', '14px')
    .attr('text-anchor', 'start')
  
  reg_builder_group.append('rect')
    .attr('transform', `translate(0,${geno_offset.y-27+row_scale(1)})`)
    .attr('x',geno_x_scale(0))
    .attr('width',geno_x_scale(num_g1)-geno_x_scale(0))
    .attr('height',row_scale.bandwidth())
    .attr("fill", "green")
    .attr("opacity",0.6)
  reg_builder_group.append('text')
    .text(num_g1)
    .attr('y',12)
    .attr('transform', `translate(${geno_x_scale(num_g1)+5},${geno_offset.y-27+row_scale(1)})`)
    .attr('font-size', '14px')
    .attr('text-anchor', 'start')
  
  // Reg Type Buttons
  const button_width = geno_panel_width - 2*(geno_offset.x)
  let rtype_group = reg_builder_group.append('g')
    .attr('transform',`translate(${geno_offset.x},${geno_offset.y+row_scale.range()[1]-10})`)
    .selectAll("empty").data(reg_types).enter()
  let button_scale = d3.scaleBand()
    .domain(d3.set(reg_types,d=>d.rtype).values())
    .range([0,75])
    .paddingInner(0.2)
  let button_color = ["grey","yellow"]
  
  rtype_group.append('rect')
    .attr('class','rtype_button')
    .attr('transform', d=> `translate(0,${button_scale(d.rtype)})`)
    .attr('width',button_width)
    .attr('height',button_scale.bandwidth())
    .attr('stroke',d3.hcl(0,0,70))
    .attr('stroke-width', 3)
    .attr("fill", d => button_color[d.s])
    .attr("opacity",0.5)
    .on('click', function(d) {
      let new_select = [0,0,0]
      new_select[d.rtype] = 1
      $0.value = z.addCol("s",new_select,reg_types)
    })
  rtype_group.append('text')
    .text(d => d.rname)
    .attr('transform', d=> `translate(${button_width/2},${button_scale(d.rtype)+15})`)
    .attr('font-size', '12px')
    .attr('text-anchor', 'middle')
    .attr('font-weight','bold')
  
  
  let cov_select_group = reg_builder_group.append('g')
    .attr('transform',`translate(${geno_panel_width},0)`)
  
  cov_selection(cov_select_group,reg_panel_width)

}
)}

function _cov_selection(z,cov_select,reg_builder_size,d3,gen_display_string,run_status,reg_types,gen_model_string,target_var,reg_list,save_covariate_data,$0,$1,$2,$3,$4,canonical,$5,hist_data){return(
function cov_selection(cov_select_group,reg_panel_width){
  
  let selected = z.getCol("vname", z.filter(r => r.s === 1, cov_select))
  // console.log(selected)
  
  // Section Outline
  cov_select_group.append('rect')
    .attr('width', reg_panel_width)
    .attr('height', reg_builder_size.height)
    .attr('stroke', d3.hcl(0,0,70))
    .attr('stroke-width', 1.5)
    .attr('fill','none')
    .lower()
  
  cov_select_group.append('text')
    .text(gen_display_string())
    .attr('text-anchor', 'middle')
    .attr('x', reg_panel_width/2)
    .attr('y', 30)
    .attr('font-size', '24px')
    .attr('font-style', 'italic')
    .attr('font-weight','bold')
    .attr('font-family','serif')
  
  // Run Button
  const button_width = 70
  const button_height = 20
  let run_button_color = ["grey","red"]
  cov_select_group.append('rect')
    .attr('class','run_button')
    .attr('transform', d=> `translate(${reg_panel_width-90},10)`)
    .attr('width',button_width)
    .attr('height',button_height)
    .attr('stroke',d3.hcl(0,0,70))
    .attr('stroke-width', 3)
    .attr("fill", run_button_color[run_status])
    .attr("opacity",0.2)
    .on('click', function(d) {
      if (run_status === 0){
        let rtype = z.filter(r => r.s ===1, reg_types)
        let model_string = gen_model_string(target_var, z.getCol("vname", z.filter(r => r.s === 1, cov_select)));
        let next_ix = reg_list.length;
        let reg_list_copy = [...reg_list]
        reg_list_copy.unshift({ cmd:"run_reg", model_str:model_string, rtype:rtype[0].rtype, save_cov:save_covariate_data})
        $0.value = [...reg_list_copy];
        $1.value = [{ cmd:"run_reg", model_str:model_string, rtype:rtype[0].rtype, save_cov:save_covariate_data}];
        $2.value = 0;
        $3.value = 1;
      }
    })
  cov_select_group.append('text')
    .text("Run")
    .attr('transform', d=> `translate(${reg_panel_width-55},25)`)
    .attr('font-size', '16px')
    .attr('text-anchor', 'middle')
    .attr('font-weight','bold')


  // Toggles for saving covariate regression info & Flipping reg equation
  let button_size = 11;
  
  let cov_button_color = 'grey';
  if (save_covariate_data === 1){
    cov_button_color = 'yellow';
  }
  cov_select_group.append('rect')
    .attr('id','Cov_Save_button')
    .attr('transform', d=> `translate(${button_size},5)`)
    .attr('width', button_size)
    .attr('height', button_size)
    .attr('stroke',d3.hcl(0,0,70))
    .attr('stroke-width', 1)
    .attr("fill", cov_button_color)
    .attr("opacity", 0.5)
    .on('click', function(d) {
      $4.value = (save_covariate_data + 1) % 2;
    })
  cov_select_group.append('text')
    .text('Save Cov Data')
    .attr('transform', d=> `translate(${(button_size)*2+5},${button_size+3})`)
    .attr('font-size', '11px')
    .attr('text-anchor', 'start')
    .attr('font-weight','bold')

  let flip_button_color = 'grey';
  if (!canonical){
    flip_button_color = 'yellow';
  }
  cov_select_group.append('rect')
    .attr('id','Flip_button')
    .attr('transform', d=> `translate(${button_size},${button_size+8})`)
    .attr('width', button_size)
    .attr('height', button_size)
    .attr('stroke',d3.hcl(0,0,70))
    .attr('stroke-width', 1)
    .attr("fill", flip_button_color)
    .attr("opacity", 0.5)
    .on('click', function(d) { $5.value = !canonical })
  cov_select_group.append('text')
    .text('Flip Equation')
    .attr('transform', d=> `translate(${(button_size)*2+5},${button_size*2 + 5})`)
    .attr('font-size', '11px')
    .attr('text-anchor', 'start')
    .attr('font-weight','bold')

  
  // Show selected covariates
  const hist_height = 75
  
  let hgroups = []
  let hgroup_scale = d3.scaleBand().domain(selected).range([0,reg_panel_width])
  const hist_width = hgroup_scale.bandwidth()*0.9;
  const hist_margin = (hgroup_scale.bandwidth() - hist_width)/2
  
  for (var i=0; i < selected.length; i++){
    let v = selected[i]
    let v_data = z.filter(d => d.var_name === v, hist_data)
    let v_data0 = z.filter(d => d.response === 0, v_data)
    let v_data1 = z.filter(d => d.response === 1, v_data)
    
    hgroups[i] = cov_select_group.append('g')
      .attr("transform", `translate(${hgroup_scale(v)},${reg_builder_size.height/3})`)
      .attr('class',v)
    
    hgroups[i].append('text')
      .text(v)
      .attr('class', v)
      .attr('x',hgroup_scale.bandwidth()/2)
      .attr('y',hist_height+20)
      .attr('font-size', '12px')
      .attr('text-anchor','middle')

    var x_scale = d3.scaleLinear()
        .domain([d3.min(v_data, d => d.xmin),d3.max(v_data, d => d.xmax)])
        .range([0, hist_width]);

    var y_scale = d3.scaleLinear()
      .domain([0, d3.max(v_data, d => d.count)]).nice()
      .range([hist_height,0])
    
    hgroups[i].append("g")
      .attr('class', v)
      .attr("transform", `translate(${hist_margin},0)`)
      .attr("fill", "purple")
      .attr("opacity",0.6)
      .selectAll("empty").data(v_data0).enter().append("rect")
        .attr("x", d => x_scale(d.xmin))
        .attr("y", d=> y_scale(d.count))
        .attr("width", d => x_scale(d.xmax) - x_scale(d.xmin))
        .attr("height",d=> hist_height - y_scale(d.count))
    
    hgroups[i].append("g")
      .attr('class', v)
      .attr("transform", `translate(${hist_margin},0)`)
      .attr("fill", "green")
      .attr("opacity",0.4)
      .selectAll("empty").data(v_data1).enter().append("rect")
        .attr("x", d => x_scale(d.xmin))
        .attr("y", d=> y_scale(d.count))
        .attr("width", d => x_scale(d.xmax) - x_scale(d.xmin))
        .attr("height",d=> hist_height - y_scale(d.count))
    
    hgroups[i].append('line')
        .attr('x1',hist_margin)
        .attr('x2',hist_margin+hist_width)
        .attr('y1',hist_height)
        .attr('y2',hist_height)
        .attr('stroke', "black")
        .attr('stroke-width',1)
  }
}
)}

function _group_legend(var_dist2_size,legend_height,d3,target_var,var_stat_scale){return(
function group_legend(legend_group){
  
  legend_group.append('rect')
    .attr('width', var_dist2_size.width)
    .attr('height', legend_height)
    .attr('stroke', d3.hcl(0,0,70))
    .attr('stroke-width', 2)
    .attr('fill','none')
    .lower()
  
  let legend_scale = d3.scaleBand()
    .domain(['genotype','pvals'])
    .range([0, var_dist2_size.width])
    .paddingInner(0.2)
    .paddingOuter(0.1)
    .align(0.5)

  const box_height = legend_height * 0.45;
  const box_offset = legend_height * 0.45;
  const txt_offset = box_offset + 12; // 12 point font
  const hdr_offset = legend_height * 0.3;
  
  // genotype colors
  let geno_data = [{name:'0',color:"purple",opacity:0.6},
                   {name:'1',color:"green",opacity:0.4}]

  let geno_group = legend_group.append('g')
    .attr('transform',`translate(${legend_scale('genotype')}, 0)`)
    .selectAll("empty").data(geno_data).enter()
  
  let geno_scale = d3.scaleBand()
    .domain(d3.set(geno_data, d=> d.name).values())
    .range([0,legend_scale.bandwidth()])
    .paddingInner(0.1)
    .align(0.5)

  geno_group.append('text')
    .text(target_var)
    .attr('font-size','12px')
    .attr('font-weight', 'bold')
    .attr("text-anchor", "middle")
    .attr('x', legend_scale.bandwidth()/2)
    .attr('y', hdr_offset)
  
  geno_group.append('rect')
    .attr('fill', d=> d.color)
    .attr('opacity',d=> d.opacity)
    .attr('width', geno_scale.bandwidth())
    .attr('height', box_height)
    .attr('x',d => geno_scale(d.name))
    .attr('y', box_offset)
  
  geno_group.append('text')
    .text(d => d.name)
    .attr('font-size','12px')
    .attr("text-anchor", "middle")
    .attr('x',d => geno_scale(d.name)+(geno_scale.bandwidth()*0.5))
    .attr('y', txt_offset)
  
  // p-value significance colors
  let pval_data = [{name:'p >= 0.05', color:var_stat_scale(0.1)},
                   {name:'p < 0.05',  color:var_stat_scale(0.01)}]

  let pval_group = legend_group.append('g')
    .attr('transform',`translate(${legend_scale('pvals')}, 0)`)
    .selectAll("empty").data(pval_data).enter()
  
  let pval_scale = d3.scaleBand()
    .domain(d3.set(pval_data, d=> d.name).values())
    .range([0,legend_scale.bandwidth()])
    .paddingInner(0.1)
    .align(0.5)

  pval_group.append('text')
    .text('Significance')
    .attr('font-size','12px')
    .attr('font-weight', 'bold')
    .attr("text-anchor", "middle")
    .attr('x', legend_scale.bandwidth()/2)
    .attr('y', hdr_offset)
  
  pval_group.append('rect')
    .attr('fill', d=> d.color)
    .attr('opacity',d=> d.opacity)
    .attr('width', pval_scale.bandwidth())
    .attr('height', box_height)
    .attr('x',d => pval_scale(d.name))
    .attr('y', box_offset)
  
  pval_group.append('text')
    .text(d => d.name)
    .attr('font-size','12px')
    .attr("text-anchor", "middle")
    .attr('x',d => pval_scale(d.name)+(pval_scale.bandwidth()*0.5))
    .attr('y', txt_offset)
}
)}

function _gen_display_string(canonical,target_var){return(
function gen_display_string(){
  var disp_str;
  if (canonical) { disp_str = "PheCode ~ " + target_var + " +" }
  else { disp_str = target_var + " ~ PheCode +" }
  return disp_str
}
)}

function _gen_model_string(canonical){return(
function gen_model_string(target, cov_list) {
  var model_string;
  if (canonical) { model_string = "PheCode~" + target }
  else { model_string = target + "~PheCode" }
  
  for (var i=0; i < cov_list.length; i++){
    model_string = model_string + "+" + cov_list[i];
  }
  
  return model_string
}
)}

function _15(md){return(
md `### Regression Panel Functions`
)}

function _reg_panel(draw_cat_legend,draw_volcano_legend,thresh_select_panel,plot_logOdds,plot_volcano,drawTable,$0){return(
function reg_panel(svg_elem) {
  
  let main_group = svg_elem.append('g')
    .attr('id', 'main')
  
  draw_cat_legend(main_group)
  draw_volcano_legend(main_group)
  thresh_select_panel(main_group)
  plot_logOdds(main_group)
  plot_volcano(main_group)
  drawTable(main_group)
  
  $0.value = 0
}
)}

function _plot_logOdds(z,thresh_value,reg_data,d3,reg_plot_size,low_lim_col,up_lim_col,reg_plot_legend_height,buffer,phe_cat_colors,vol_color,xaxis_name,phe_cat_markers,point_col,model_type,$0){return(
function plot_logOdds(main_group)  {
  // NOTE 
  // This is is not necessarily a log odds plot anymore!! There is a toggle to switch between Odds Ratio and Log Odds Ratio.
  // Toggle button code is at the bottom of this cell

  let margin = 15
  let amargin = 50
  // filter data
  let filtered_reg = z.filter(r => r.pval < thresh_value, reg_data)
  // define scales
  let logodds_y_scale = d3.scaleBand()
                          .domain(d3.set(filtered_reg, d => d.PheCode).values())
                          .range([margin,reg_plot_size.height-amargin])
  let logodds_x_scale = d3.scaleLinear()
                          .domain([d3.min(filtered_reg, d=>d[low_lim_col]), 
                                   d3.max(filtered_reg, d=>d[up_lim_col])])
                          .range([margin,reg_plot_size.width-margin])
                          .nice()
  
  // plot groups
  let log_odds_group = main_group.append('g')
    .attr('transform', 'translate(0,'+reg_plot_legend_height+')')
    .attr('id', 'log_odd_plot') 
  
  let beta_group = log_odds_group.selectAll('empty')
    .attr('transform', `translate(${buffer},${buffer})`)
    .data(filtered_reg)
    .enter()
  
  // Background
  log_odds_group.append('rect')
    .attr('width', reg_plot_size.width)
    .attr('height', reg_plot_size.height)
    .attr('stroke', d3.hcl(0,0,70))
    .attr('stroke-width', 2)
    .attr('fill','transparent')
    .lower()
    .on('click', function(d) {
      main_group.selectAll('.logodds')
        .attr('opacity', 1.0)
        .attr('stroke', d => phe_cat_colors(d.category_id))
      main_group.selectAll('.volcano')
        .attr('opacity', 1.0)
        .attr('stroke', d => vol_color(d.pval))
      main_group.selectAll('.phelabel_click').remove()
      main_group.selectAll('.phecode_row').attr("fill","grey")
    })
  
  if (reg_data.length == 1 && reg_data[0].msg === "no_data"){
    log_odds_group.append('text')
      .attr('class','placeholder')
      .text('Effect Size Plot')
      .attr('text-anchor','middle')
      .attr('transform',`translate(${reg_plot_size.width/2},${reg_plot_size.height/3})`)
      .attr('font-size','18px')
      .attr('font-weight','bold')
  }
  else {
    // axes
    log_odds_group.append('line')
      .attr('x1', d=>logodds_x_scale(0))
      .attr('y1', d=>logodds_y_scale.range()[0])
      .attr('x2', d=>logodds_x_scale(0))
      .attr('y2', d=>logodds_y_scale.range()[1])
      .attr('stroke',d => "black")
      .attr('stroke-width', 1)

    let log_xaxis = d3.axisBottom(logodds_x_scale).ticks(6)

    log_odds_group.append('g')
      .attr('id', 'xaxis')
      .attr('transform', `translate(0,${logodds_y_scale.range()[1]})`)
      .call(log_xaxis)

    log_odds_group.append('text')
      .text(xaxis_name)
      .attr('id', 'xaxis_label')
      .attr('transform', `translate(${reg_plot_size.width/2},${logodds_y_scale.range()[1]+35})`)
      .attr('text-anchor', 'middle')
      .attr('font-size', '13px')
      .attr('font-weight', 'bold')

    // confidence intervals
    beta_group.append('line')
      .attr('class','logodds')
      .attr('id',d=>d.Pheno_id)
      .attr('x1', d=>logodds_x_scale(d[low_lim_col]))
      .attr('y1', d=>logodds_y_scale(d.PheCode))
      .attr('x2', d=>logodds_x_scale(d[up_lim_col]))
      .attr('y2', d=>logodds_y_scale(d.PheCode))
      .attr('stroke',d => phe_cat_colors(d.category_id))
      .attr('stroke-width', 2)

    // point estimate 
    beta_group.append('path')
      .attr('id',d=>d.Pheno_id)
      .attr('class','logodds')
      .attr('fill', d => phe_cat_colors(d.category_id))
      .attr("d", d => phe_cat_markers(d.category_id)) 
      .attr("transform",  d => `translate(${logodds_x_scale(d[point_col])}, ${logodds_y_scale(d.PheCode)})`)
      .attr('x', d=>logodds_x_scale(d[point_col])) // path doesn't use x/y
      .attr('y', d=>logodds_y_scale(d.PheCode)) // but this makes it easy to place text on click
      .on('mouseover', function(d) {
        log_odds_group.append('text')
          .attr("class","phelabel")
          .text(d.Phenotype)
          .attr('text-anchor', 'start')
          .attr('x', +d3.select(this).attr('x')+10)
          .attr('y', +d3.select(this).attr('y')-10)
          .attr('font-weight', 'bold')
          .style('font-size', 13)
          .raise()
      })
      .on('mouseout', (d, i, arr) => {
        log_odds_group.selectAll('.phelabel').remove()
      })
      .on('click', function(d) {
        // reset attributes
        main_group.selectAll('.logodds')
          .attr('opacity', 0.2)
          .attr('stroke', d => phe_cat_colors(d.category_id))
        main_group.selectAll('.volcano')
          .attr('opacity', 0.2)
          .attr('stroke', d => vol_color(d.pval))
        main_group.selectAll('.phelabel_click').remove()
        main_group.selectAll('.phecode_row').attr("fill","grey")

        console.log(d.Pheno_id)
        // highlight selected point in both plots
        log_odds_group.selectAll('#'+d.Pheno_id)
          .attr('opacity', 1.0)
          .attr('stroke','black')
        main_group.selectAll('.volcano')
          .filter(function(r) { 
            return (r.Pheno_id === d.Pheno_id)
          })
          .attr('opacity', 1.0)
          .attr('stroke','black')
        // draw text 
        log_odds_group.append('text')
          .attr("class","phelabel_click")
          .text(d.Phenotype)
          .attr('x',+d3.select(this).attr('x')+10)
          .attr('y', +d3.select(this).attr('y')-10)
          .attr('font-weight', 'bold')
          .style('font-size', 13)
          .attr('text-anchor', 'start')
          .raise()
        // highlight selected point in table (if it's shown)
        // color row
        main_group.selectAll('.phecode_row')
          .filter(function(r) { 
            return (r.Pheno_id === d.Pheno_id)
          })
          .attr("fill",d => phe_cat_colors(d.category_id))
      })

    
  if (xaxis_name !== "Regression Coefficient"){
    // Toggle for Odds Ratio display
    let button_size = 15;
    let button_color = 'grey';
    if (xaxis_name === 'Odds Ratio'){
      button_color = 'yellow';
    }
      
    log_odds_group.append('rect')
      .attr('id','OR_button')
      .attr('transform', d=> `translate(${margin},${reg_plot_size.height-button_size-5})`)
      .attr('width', button_size)
      .attr('height', button_size)
      .attr('stroke',d3.hcl(0,0,70))
      .attr('stroke-width', 3)
      .attr("fill", button_color)
      .attr("opacity", 0.5)
      .on('click', function(d) {
        if (model_type === 0){
          $0.value = 1;
        }
        else{
          $0.value = 0;
        }
      })
    log_odds_group.append('text')
      .text('Toggle OR')
      .attr('transform', d=> `translate(${margin+button_size+5},${reg_plot_size.height-button_size+7})`)
      .attr('font-size', '11px')
      .attr('text-anchor', 'start')
      .attr('font-weight','bold')
  }

    return beta_group
 }
}
)}

function _plot_volcano(d3,reg_data,reg_plot_size,point_col,reg_plot_legend_height,buffer,phe_cat_colors,vol_color,xaxis_name){return(
function plot_volcano(main_group)  {
  let margin = 15
  let amargin = 50
  // define scales
  let volcano_y_scale = d3.scaleLinear()
                          .domain(d3.extent(reg_data, d=>d.neg_log_p))
                          .range([reg_plot_size.height-amargin,margin])
                          .nice()
  let volcano_x_scale = d3.scaleLinear()
                          .domain(d3.extent(reg_data, d=>d[point_col]))
                          .range([amargin,reg_plot_size.width-margin])
                          .nice()
  
  // plot groups
  let volcano_group = main_group.append('g')
    .attr('transform', 
          'translate('+(reg_plot_size.width)+','+reg_plot_legend_height+')')
    .attr('id', 'volcano_plot')
  
  let beta_group = volcano_group.selectAll('empty')
    .attr('transform', `translate(${buffer},${buffer})`)
    .data(reg_data)
    .enter()
  
  // Background
  volcano_group.append('rect')
    .attr('width', reg_plot_size.width)
    .attr('height', reg_plot_size.height)
    .attr('stroke', d3.hcl(0,0,70))
    .attr('stroke-width', 2)
    .attr('fill','transparent')
    .lower()
    .on('click', function(d) {
      main_group.selectAll('.logodds')
        .attr('opacity', 1.0)
        .attr('stroke', d => phe_cat_colors(d.category_id))
      main_group.selectAll('.volcano')
        .attr('opacity', 1.0)
        .attr('stroke', d => vol_color(d.pval))
      main_group.selectAll('.phelabel_click').remove()
      main_group.selectAll('.phecode_row').attr("fill","grey")
    })
  
  if (reg_data.length == 1 && reg_data[0].msg === "no_data"){
    volcano_group.append('text')
      .attr('class','placeholder')
      .text('Volcano Plot')
      .attr('text-anchor','middle')
      .attr('transform',`translate(${reg_plot_size.width/2},${reg_plot_size.height/3})`)
      .attr('font-size','18px')
      .attr('font-weight','bold')
  }
  else {
    // axes
    let vol_xaxis = d3.axisBottom(volcano_x_scale).ticks(6)
    let vol_yaxis = d3.axisLeft(volcano_y_scale).ticks(6)

    volcano_group.append('g')
      .attr('id', 'yaxis')
      .attr('transform', 'translate(0,'+volcano_y_scale.range()[0]+')')
      .call(vol_xaxis)
    volcano_group.append('text')
      .text('-log(p)')
      .attr('id', 'yaxis_label')
      .attr('transform', `rotate(270) translate(-${reg_plot_size.height/2},20)`)
      .attr('text-anchor', 'middle')
      .attr('font-size', '15px')
      .attr('font-weight', 'bold')

    volcano_group.append('g')
      .attr('id', 'xaxis')
      .attr('transform', 'translate('+(amargin)+',0)')
      .call(vol_yaxis)
    volcano_group.append('text')
      .text(xaxis_name)
      .attr('id', 'xaxis_label')
      .attr('transform', `translate(${reg_plot_size.width/2},${volcano_y_scale.range()[0]+35})`)
      .attr('text-anchor', 'middle')
      .attr('font-size', '13px')
      .attr('font-weight', 'bold')

    // point estimate 
    beta_group.append('circle')
      .attr('id',d=>d.Pheno_id)
      .attr("class","volcano")
      .attr('cx', d=>volcano_x_scale(d[point_col]))
      .attr('cy', d=>volcano_y_scale(d.neg_log_p))
      .attr('r', 3)
      .attr('stroke-width', 2)
      .attr('stroke', d => vol_color(d.pval))
      .attr('fill', d => vol_color(d.pval))
      .on('mouseover', function(d) {
        const cx = +d3.select(this).attr('cx')
        const phe_len = d3.select(this).data()[0].Phenotype.length
        let ta = 'start'
        if (cx+(phe_len*6) > volcano_x_scale.range()[1]){
          ta = 'end'
        }
        volcano_group.append('text')
          .attr("class","phelabel")
          .text(d.Phenotype)
          .attr('x',cx+5)
          .attr('y', +d3.select(this).attr('cy')-10)
          .attr('font-weight', 'bold')
          .style('font-size', 13)
          .attr('text-anchor', ta)
          .raise()
      })
      .on('mouseout', () => {
        volcano_group.selectAll('.phelabel').remove()
      })
      .on('click', function(d) {
        // reset attributes
        main_group.selectAll('.logodds')
          .attr('opacity', 0.2)
          .attr('stroke', d => phe_cat_colors(d.category_id))
        main_group.selectAll('.volcano')
          .attr('opacity', 0.2)
          .attr('stroke', d => vol_color(d.pval))
        main_group.selectAll('.phelabel_click').remove()
        main_group.selectAll('.phecode_row').attr("fill","grey")

        // highlight selected point in both plots
        volcano_group.selectAll('#'+d.Pheno_id)
          .attr('opacity', 1.0)
          .attr('stroke','black')
        main_group.selectAll('.logodds')
          .filter(function(r) { 
            return (r.Pheno_id === d.Pheno_id)
          })
          .attr('opacity', 1.0)
          .attr('stroke','black')
        // draw text 
        const cx = +d3.select(this).attr('cx')
        const phe_len = d3.select(this).data()[0].Phenotype.length
        let ta = 'start'
        if (cx+(phe_len*6) > volcano_x_scale.range()[1]){
          ta = 'end'
        }
        volcano_group.append('text')
          .attr("class","phelabel_click")
          .text(d.Phenotype)
          .attr('x',cx+5)
          .attr('y', +d3.select(this).attr('cy')-10)
          .attr('font-weight', 'bold')
          .style('font-size', 13)
          .attr('text-anchor', ta)
          .raise()
        // highlight selected point in table (if it's shown)
        // color row
        main_group.selectAll('.phecode_row')
          .filter(function(r) { 
            return (r.Pheno_id === d.Pheno_id)
          })
          .attr("fill",d => phe_cat_colors(d.category_id))
      })
  }
}
)}

function _draw_volcano_legend(reg_plot_legend_height,cat_legend_width,d3){return(
function draw_volcano_legend(main_group){
  const legend_width = 90
  const legend_height = reg_plot_legend_height
  const entries = [{type:'Bonferroni',color:"gold"},
                   {type:'FDR',color:"blue"},
                   {type:'Insignificant',color:"gray"}]
  
  let vol_legend = main_group.append('g')
    .attr('transform',`translate(${cat_legend_width},0)`)
  
  // outline
  vol_legend.append('rect')
    .attr('width', legend_width)
    .attr('height', legend_height)
    .attr('stroke', d3.hcl(0,0,70))
    .attr('stroke-width', 2)
    .attr('fill','white')
    .lower()
  
  // scale
  let text_scale = d3.scaleBand()
    .domain(d3.set(entries, d=> d.type).values())
    .range([0,legend_height-25])
    .paddingInner(0.05)
  
    // Title
  vol_legend.append('text').text('Thresholds')
    .attr('text-anchor', 'middle')
    .attr('x', (legend_width/2))
    .attr('y', 15)
    .attr('font-size', '13px')
    .attr('font-weight', 'bold')
  
  let legend_data_group = vol_legend.selectAll('g').data(entries).enter().append('g')
  
  legend_data_group.append('circle')
    .attr('cx',10)
    .attr('cy',d =>  30 + text_scale(d.type)).attr('r',3)
    .attr('fill', d => d.color)
    .attr('stroke-width', 2)
    .attr('stroke', d => d.color)
  
  legend_data_group.append("text").text(d => d.type)
		.attr('x', 20)
    .attr('y', d =>  35 + text_scale(d.type))
		.attr('text-anchor', 'start')
    .attr('font-size', '12px')
  
}
)}

function _draw_cat_legend(cat_legend_width,reg_plot_legend_height,d3,phe_categories,phe_cat_colors,phe_cat_markers){return(
function draw_cat_legend(main_group){
  
  let cat_square_size = 5
  let margin = 10
  let squares_per_row = 5
  let y_space = 25
  
  let legend_group = main_group.append('g')
    .attr('id', 'reg_cat_legend')
    .attr('transform',`translate(${margin},${margin})`)
  
  //outline
  main_group.append('rect')
    .attr('width', cat_legend_width)
    .attr('height', reg_plot_legend_height)
    .attr('stroke', d3.hcl(0,0,70))
    .attr('stroke-width', 2)
    .attr('fill','none')
    .lower()
  // description
  main_group.append('text')
    .text("PheCode Categories")
    .attr('text-anchor', 'middle')
    .attr('x', cat_legend_width*0.75)
    .attr('y', reg_plot_legend_height*0.9)
    .attr('font-size', '15px')
    .attr('font-weight', 'bold')
  
  // category legend
  let square_scale = d3.scaleBand()
    .domain(d3.range(squares_per_row))
    .range([0,cat_legend_width-20])
    .paddingInner(0.1)

  let annotation_group = legend_group.selectAll("empty").data(phe_categories).enter()

  annotation_group.append("path")
    .attr('transform', (d,i) => {
      var y_off = Math.floor(i/squares_per_row)
      var square_offset_y = y_off*y_space
      var x_off = i % squares_per_row
      var square_offset_x = square_scale(x_off)
      return 'translate('+square_offset_x+','+square_offset_y+')'
    })
   .attr('fill', d => phe_cat_colors(d.category))
   .attr("d", r => phe_cat_markers(r.category)) 
 
  
  annotation_group.append("text")
    .attr('transform', (d,i) => {
      var y_off = Math.floor(i/squares_per_row)
      var square_offset_y = y_off*y_space
      var x_off = i % squares_per_row
      var square_offset_x = square_scale(x_off)
      return 'translate('+square_offset_x+','+square_offset_y+')'
    })
    .text(d => d.category_string)
    .attr('x', cat_square_size+3)
    .attr('y', cat_square_size)
    .style('font-size', 11)
    .attr('text-anchor', 'start')

}
)}

function _thresh_select_panel(reg_panel_size,cat_legend_width,reg_plot_legend_height,d3,thresh_types,$0,z){return(
function thresh_select_panel(main_group){
  
  let ts_width = reg_panel_size - (cat_legend_width+90)
  let mw = 10
  
  let thresh_group = main_group.append('g')
    .attr('transform',`translate(${(cat_legend_width+90)},0)`)
  
  thresh_group.append('rect')
    .attr('width', ts_width)
    .attr('height', reg_plot_legend_height)
    .attr('stroke', d3.hcl(0,0,70))
    .attr('stroke-width', 2)
    .attr('fill','none')
    .lower()
  
  // Description
  thresh_group.append('text').text('Threshold Selection')
    .attr('text-anchor', 'middle')
    .attr('x', (ts_width/2))
    .attr('y', 15)
    .attr('font-size', '13px')
    .attr('font-weight', 'bold')
  
  // Thresh Type Buttons
  let tb_group = thresh_group.append('g')
    .selectAll("empty").data(thresh_types).enter()
  
  const button_width = ts_width - 2*(mw)
  let button_scale = d3.scaleBand()
    .domain(d3.set(thresh_types, d=> d.name).values())
    .range([0,reg_plot_legend_height-40])
    .paddingInner(0.2)
  let button_color = ["grey","yellow"]
  
  tb_group.append('rect')
    .attr('class','ts_button')
    .attr('transform', d=> `translate(${mw},${button_scale(d.name)+25})`)
    .attr('width',button_width)
    .attr('height',button_scale.bandwidth())
    .attr('stroke',d3.hcl(0,0,70))
    .attr('stroke-width', 3)
    .attr("fill", d => button_color[d.s])
    .attr("opacity",0.5)
    .on('click', function(d) {
      let new_select = [0,0]
      new_select[d.id] = 1
      $0.value = z.addCol("s",new_select,thresh_types)
    })
  tb_group.append('text')
    .text(d => d.name)
    .attr('transform', d=> `translate(${ts_width/2},${button_scale(d.name)+42})`)
    .attr('font-size', '12px')
    .attr('text-anchor', 'middle')
    .attr('font-weight','bold')

  
  
}
)}

function _drawTable(reg_plot_legend_height,reg_plot_size,reg_panel_size,reg_table_height,d3,reg_data,model_string,table_data,table_cols,z,point_col,phe_cat_colors,vol_color,buffer,start_ix,num_rows,$0){return(
function drawTable(main_group) {
  let table_group = main_group.append('g')
    .attr('transform', 
          'translate(0,'+(reg_plot_legend_height+reg_plot_size.height)+')')
    .attr('id', 'reg_table')
  
  table_group.append('rect')
    .attr('width', reg_panel_size)
    .attr('height', reg_table_height)
    .attr('stroke', d3.hcl(0,0,70))
    .attr('stroke-width', 2)
    .attr('fill','none')
    .lower()
  
  if (reg_data.length == 1 && reg_data[0].msg === "no_data"){
    table_group.append('text')
      .attr('class','placeholder')
      .text('Regression Data Table')
      .attr('text-anchor','middle')
      .attr('transform',`translate(${reg_panel_size/2},${reg_table_height/3})`)
      .attr('font-size','18px')
      .attr('font-weight','bold')
  }
  else {
    // show model string for current regression results
    table_group.append('text')
      .text(model_string)
      .attr('text-anchor','middle')
      .attr('transform',`translate(${reg_panel_size/2},${20})`)
      .attr('font-size','18px')
      .attr('font-style', 'italic')
      .attr('font-weight','bold')
      .attr('font-family','serif')

    let mi = 40;
    let table_width = 0.85*reg_panel_size
    let f_pval = d3.format(".3e") // format pvalues
    let f_beta = d3.format(".4f") // format beta values

    let rows = d3.set(table_data,d=>d.PheCode).values()
    rows.unshift('header')
    let row_scale = d3.scaleBand()
      .domain(rows)
      .range([mi,reg_table_height])
      .paddingInner(0.1)
      .paddingOuter(0.05)

    // header row
    let header_group = table_group.append('g')
      .attr('transform',`translate(${mi},${row_scale('header')})`)
      .selectAll("empty").data(table_cols).enter()

    header_group.append('text')
      .text(d => d.name)
      .attr('font-size','12px')
      .attr('font-weight','bold')
      .attr('x',d=> d.offset)

    // data rows
    let data_group = table_group.append('g')
      .selectAll("empty").data(table_data).enter()
      .append('g').attr('transform', d => `translate(${mi},${row_scale(d.PheCode)})`)

    // PheCode
    let col = z.filter(r => r.name === 'PheCode',table_cols)[0]
    data_group.append('text')
      .attr('class','phecode_row_text')
      .attr('id',d => d.Pheno_id)
      .text(d => d.PheCode)
      .attr('font-size','12px')
      .attr('x', col.offset+(col.width)*0.8)
      .attr('text-anchor','end')

    // Phenotype
    col = z.filter(r => r.name === 'Phenotype',table_cols)[0]
    data_group.append('text')
      .attr('class','phecode_row_text')
      .attr('id',d => d.Pheno_id)
      .text(d => d.Phenotype)
      .attr('font-size','12px')
      .attr('x', col.offset)

    // Count
    col = z.filter(r => r.name === 'Count',table_cols)[0]
    data_group.append('text')
      .attr('class','phecode_row_text')
      .attr('id',d => d.Pheno_id)
      .text(d => d.count)
      .attr('font-size','12px')
      .attr('x', col.offset)

    // OR
    col = z.filter(r => r.name === point_col,table_cols)[0]
    data_group.append('text')
      .attr('class','phecode_row_text')
      .attr('id',d => d.Pheno_id)
      .text(d => f_beta(d[point_col]))
      .attr('font-size','12px')
      .attr('x', col.offset)

    // p-value
    col = z.filter(r => r.name === 'P-value',table_cols)[0]
    data_group.append('text')
      .attr('class','phecode_row_text')
      .attr('id',d => d.Pheno_id)
      .text(d => f_pval(d.pval))
      .attr('font-size','12px')
      .attr('x', col.offset)

    // category
    col = z.filter(r => r.name === 'Category',table_cols)[0]
    data_group.append('text')
      .attr('class','phecode_row_text')
      .attr('id',d => d.Pheno_id)
      .text(d => d.category)
      .attr('font-size','12px')
      .attr('x', col.offset)

    // row background
    data_group.append('rect')
      .attr('class','phecode_row')
      .attr('id',d => d.Pheno_id)
      .attr('y',-row_scale.bandwidth()+3)
      .attr('height',row_scale.bandwidth())
      .attr('width',col.offset)
      .attr('fill',"grey")
      .attr("opacity",0.2)
      .on("click", function(d){
        // reset
        main_group.selectAll('.logodds')
          .attr('opacity', 0.2)
          .attr('stroke', d => phe_cat_colors(d.category_id))
        main_group.selectAll('.volcano')
          .attr('opacity', 0.2)
          .attr('stroke', d => vol_color(d.pval))
        main_group.selectAll('.phelabel_click').remove()
        data_group.selectAll('.phecode_row').attr("fill","grey")
        // highlight selected point in both plots
        main_group.selectAll('.volcano')
          .filter(function(r) { 
            return (r.Pheno_id === d.Pheno_id)
          })
          .attr('opacity', 1.0)
          .attr('stroke','black')
        main_group.selectAll('.logodds')
          .filter(function(r) { 
            return (r.Pheno_id === d.Pheno_id)
          })
          .attr('opacity', 1.0)
          .attr('stroke','black')
        // color row
        data_group.selectAll('.phecode_row')
          .filter(function(r) { 
            return (r.Pheno_id === d.Pheno_id)
          })
          .attr("fill",d => phe_cat_colors(d.category_id))
      })
    data_group.append('rect')
      .attr('x',col.offset)
      .attr('y',-row_scale.bandwidth()+3)
      .attr('height',row_scale.bandwidth())
      .attr('width',col.width)
      .attr('fill', d=> phe_cat_colors(d.category_id))
      .attr("opacity",0.2)
      .on("click", function(d){
        // reset
        main_group.selectAll('.logodds')
          .attr('opacity', 0.2)
          .attr('stroke', d => phe_cat_colors(d.category_id))
        main_group.selectAll('.volcano')
          .attr('opacity', 0.2)
          .attr('stroke', d => vol_color(d.pval))
        main_group.selectAll('.phelabel_click').remove()
        data_group.selectAll('.phecode_row').attr("fill","grey")
        // highlight selected point in both plots
        main_group.selectAll('.volcano')
          .filter(function(r) { 
            return (r.Pheno_id === d.Pheno_id)
          })
          .attr('opacity', 1.0)
          .attr('stroke','black')
        main_group.selectAll('.logodds')
          .filter(function(r) { 
            return (r.Pheno_id === d.Pheno_id)
          })
          .attr('opacity', 1.0)
          .attr('stroke','black')
        // color row
        data_group.selectAll('.phecode_row')
          .filter(function(r) { 
            return (r.Pheno_id === d.Pheno_id)
          })
          .attr("fill",d => phe_cat_colors(d.category_id))
      })

    // scrolling arrows
    let arrow_panel_width = (reg_panel_size - (table_width*1.05))
    let arrow_width = arrow_panel_width*0.4
    let arrow_height = arrow_width*0.6
    let arrow_offset = arrow_panel_width*0.1
    let arrow_group = table_group.append('g')
      .attr('transform',`translate(${mi+table_width+buffer+arrow_offset},${reg_table_height/2 + buffer})`)

    var down_arrow = [{"x":0, "y":0},
                      {"x":arrow_width, "y":0},
                      {"x":arrow_width/2, "y":arrow_height}];
    arrow_group.selectAll("empty")
      .data([down_arrow])
      .enter().append("polygon")
        .attr('transform',`translate(0,${buffer})`)
        .attr('stroke',d3.hcl(0,0,70))
        .attr('stroke-width', 3)
        .attr("fill", "grey")
        .attr("opacity",0.2)
        .attr("points", function(d)  {
          return d.map(function(d)  {
            return [d.x, d.y].join(",");
          }).join(" ");
        })
        .on("click",function(d){
          if (!((start_ix+num_rows) >= (reg_data.length-1))) {
            $0.value = start_ix+num_rows
          }
          else{
            d3.select(this).attr("fill","black")
          }
        })

    var up_arrow = [{"x":0, "y":0},
                    {"x":arrow_width, "y":0},
                    {"x":arrow_width/2, "y":-arrow_height}];
    arrow_group.selectAll("empty")
      .data([up_arrow])
      .enter().append("polygon")
        .attr('transform',`translate(0,${-buffer})`)
        .attr('stroke',d3.hcl(0,0,70))
        .attr('stroke-width', 3)
        .attr("fill", "grey")
        .attr("opacity",0.2)
        .attr("points", function(d)  {
          return d.map(function(d)  {
            return [d.x, d.y].join(",");
          }).join(" ");
        }) 
        .on("click",function(d){
          if (start_ix >= num_rows) {
            $0.value = start_ix-num_rows
          }
          else{
            d3.select(this).attr("fill","black")
          }
        })
  }
}
)}

function _23(md){return(
md `### Threshold Calculations`
)}

function _fdr_thresh(z,reg_data,alpha)
{
	//
	// Calculate the false discovery rate threshold.
	// :param p_values: a list of p-values obtained by executing the regression
	// :param alpha: the uncorrected significance level being used (usually 0.05)
	// :type p_values: numpy array
	// :type alpha: float
	// :returns: the false discovery rate
	// :rtype: float
	///
  let p_values = z.getCol("pval",z.sortByCol("pval", "asc", reg_data));
  // console.log(z.max(p_values))
  let thresh = 0;
  const num_pvals = p_values.length;
  for(var i=0; i<num_pvals; i++){
	  let p_crit = alpha * (i+1) / num_pvals
    // console.log([i, p_crit])
		if (p_values[i] <= p_crit) {continue}
		else{
      thresh = p_values[i];
      break
    }
   }
                                
	return thresh
}


function _bon_thresh(z,reg_data,alpha)
{
  //
	// Calculate the bonferroni correction threshold.
	// Divide the power by the sum of all finite values (all non-nan values).
	// :param p_values: a list of p-values obtained by executing the regression
	// :param alpha: the uncorrected significance level being used (usually 0.05)
	// :type p_values: numpy array
	// :type alpha: float
	// :returns: The bonferroni correction
	// :rtype: float
	//
  let p_values = z.getCol("pval",z.sortByCol("pval", "asc", reg_data));
	return alpha / p_values.length
}


function _thresh_value(z,thresh_types,fdr_thresh,bon_thresh)
{
  let thresh = z.filter(r => r.s === 1, thresh_types)
  if (thresh[0].name === "FDR"){
    return fdr_thresh
  }
  else{
    return bon_thresh
  }
}


function _27(md){return(
md `### Color Scales`
)}

function _phe_cat_colors(){return(
function phe_cat_colors(category_id) {
  let remainder = category_id % 5;
  if (remainder === 0){
    return "#ff000d"
  }
  else if (remainder === 1){
    return "#f5bf03"
  }
  else if (remainder === 2){
    return "#39ad48" // "#13eac9"
  }
  else if (remainder === 3){
    return "#0343df" // "#004577"
  }
  
  return "#a00498" // "#94568c"
}
)}

function _phe_cat_markers(d3){return(
function phe_cat_markers(category_id) {
  let quotient = Math.floor(category_id / 5);
  if (quotient === 0){
    return d3.symbol().type(d3.symbolCircle).size(60)();
  }
  else if (quotient === 1){
    return d3.symbol().type(d3.symbolDiamond).size(60)();
  }
  else if (quotient === 2){
    return d3.symbol().type(d3.symbolSquare).size(75)();
  }
  
  return d3.symbol().type(d3.symbolTriangle).size(60)();
}
)}

function _var_stat_scale(d3){return(
function var_stat_scale(pval) {
  if (pval < 0.05){
    return "gold"
  }
  return d3.hcl(0,0,80)
}
)}

function _vol_color(bon_thresh,fdr_thresh){return(
function vol_color(pval) {
  if (pval < bon_thresh){
    return "gold"
  }
  if (pval < fdr_thresh){
    return "blue"
  }
  return "grey"
}
)}

function _32(md){return(
md `### Parameters`
)}

function _buffer(){return(
10
)}

function _margin(){return(
{x: 10, y: 10}
)}

function _alpha(){return(
0.05
)}

function _36(md){return(
md `#### Regression Panel`
)}

function _reg_panel_size(){return(
900
)}

function _reg_plot_size(reg_panel_size,reg_plot_legend_height,reg_table_height){return(
{width: (reg_panel_size)/2,
                  height: (reg_panel_size-reg_plot_legend_height - reg_table_height)}
)}

function _reg_plot_legend_height(){return(
100
)}

function _reg_table_height(){return(
175
)}

function _cat_legend_width(reg_panel_size){return(
reg_panel_size*(3/4)
)}

function _num_rows(){return(
7
)}

function _point_col(OR_type,model_type)
{
  return OR_type[model_type].point
}


function _up_lim_col(OR_type,model_type)
{
  return OR_type[model_type].up_lim
}


function _low_lim_col(z,OR_type,model_type)
{
  let selected_data = z.filter(r => r.s===1, OR_type)
  return OR_type[model_type].low_lim
}


function _xaxis_name(OR_type,model_type)
{
  return OR_type[model_type].name
}


function _47(md){return(
md`#### Group Panel`
)}

function _legend_height(){return(
39
)}

function _group_panel_size(){return(
700
)}

function _var_dist_size(group_panel_size){return(
{width:group_panel_size*(2/5), height:group_panel_size*(3/4)}
)}

function _var_dist2_size(group_panel_size,legend_height){return(
{width:group_panel_size*(3/5), height:group_panel_size*(3/4)-legend_height}
)}

function _reg_builder_size(group_panel_size){return(
{width:group_panel_size, height:(group_panel_size)*(1/4)}
)}

function _var1_select(z,var_comp)
{
  let num_selected = z.filter(r => r.s===1, var_comp).length
  if (num_selected > 0)
    return z.filter(r => r.s===1, var_comp)[0].vname
  else{
    return ''
  }
}


function _var2_select(z,var_comp)
{
  let num_selected = z.filter(r => r.s===1, var_comp).length
  if (num_selected === 1){
    return z.filter(r => r.s===1, var_comp)[0].vname
  }
  else if (num_selected === 2){
    return z.filter(r => r.s===1, var_comp)[1].vname
  }
  else{
    return ''
  }
}


function _55(md){return(
md `## Data`
)}

function _56(md){return(
md `### Mutables`
)}

function _target_var(response_var){return(
response_var[0].msg
)}

function _group_vars(z,group_var_init){return(
z.getCol("msg", group_var_init)
)}

function _group_data_args(){return(
[{cmd:"group_data", response: "INIT"}]
)}

function _hist_args(){return(
[{cmd:"histograms", response: "INIT"}]
)}

function _reg_args(){return(
[{ cmd:"run_reg", model_str:"", rtype:-1, save_cov:0}]
)}

function _reg_list(){return(
[{ cmd:"run_reg", model_str:"", rtype:0}]
)}

function _var_comp(group_vars)
{
  let var_info = []
  for (var i=0; i < group_vars.length; i++){
    var_info[i] = {vname: group_vars[i],
                   s:0,
                   id:i
                  }
  }
  return var_info
}


function _cov_select(group_vars)
{
  let var_info = []
  for (var i=0; i < group_vars.length; i++){
    var_info[i] = {vname: group_vars[i],
                   s:0,
                   id:i
                  }
  }
  return var_info
}


function _reg_types(){return(
[{rname:'Binary',rtype:0,s:1},
                     {rname:'Count',rtype:1,s:0},
                     {rname:'Duration',rtype:2,s:0}]
)}

function _thresh_types(){return(
[{name:'FDR',id:0,s:1},
                        {name:'Bonferroni',id:1,s:0}]
)}

function _OR_type(){return(
[{name:'Log Odds Ratio',point:'LOR',low_lim:'beta_ci_low',up_lim:'beta_ci_up'},
                   {name:'Odds Ratio',point:'OR',low_lim:'OR_ci_low',up_lim:'OR_ci_up'},
                   {name:'Regression Coefficient',point:'Beta',low_lim:'beta_ci_low',up_lim:'beta_ci_up'}]
)}

function _start_ix(){return(
0
)}

function _run_status(){return(
0
)}

function _save_covariate_data(){return(
0
)}

function _canonical(){return(
true
)}

function _model_type(){return(
0
)}

function _model_string(reg,$0)
{
  let stat_info = reg[reg.length-1]
  if (stat_info['Pheno_id'] === 'linear') {
    $0.value = 2;
  }
  else {
    $0.value = 0;
  }
  return stat_info['Phenotype']
}


function _74(md){return(
md `### Import Data`
)}

async function* _response_var($,Promises)
{
  let data_for_server = { cmd:"init_response"};

    let next_data = await $.ajax({
      url: 'http://127.0.0.1:5000/grab_data',
      dataType: 'json',
      data: JSON.stringify(data_for_server),
      contentType: 'application/json;charset=UTF-8',
      type: 'POST'
    });

    yield Promises.delay(1000, JSON.parse(next_data));
}


async function* _group_var_init($,Promises)
{
  let data_for_server = { cmd:"init_groupvars"};

    let next_data = await $.ajax({
      url: 'http://127.0.0.1:5000/grab_data',
      dataType: 'json',
      data: JSON.stringify(data_for_server),
      contentType: 'application/json;charset=UTF-8',
      type: 'POST'
    });

    yield Promises.delay(1000, JSON.parse(next_data));
}


async function* _group_data(group_data_args,$,Promises)
{
  let data_for_server = group_data_args[0];

    let next_data = await $.ajax({
      url: 'http://127.0.0.1:5000/grab_data',
      dataType: 'json',
      data: JSON.stringify(data_for_server),
      contentType: 'application/json;charset=UTF-8',
      type: 'POST'
    });

    yield Promises.delay(1000, JSON.parse(next_data));
}


async function* _reg(reg_args,$,Promises)
{
  let data_for_server = reg_args[0];

    let next_data = await $.ajax({
      url: 'http://127.0.0.1:5000/grab_data',
      dataType: 'json',
      data: JSON.stringify(data_for_server),
      contentType: 'application/json;charset=UTF-8',
      type: 'POST'
    });
  
    yield Promises.delay(1000, JSON.parse(next_data));
}


function _reg_data(z,reg)
{
  // small pvals sometimes get convereted to 0 by jsonify - this fixes them
  const pvals = z.deriveCol((r) => parseFloat(r.pval_str), reg)
  let data = z.addCol("pval", pvals, reg)
  return z.sortByCol("category_id", "asc", z.filter(r => r.PheCode != "STAT_INFO", data))
}


function _phe_categories(){return(
[{"category":0,"category_string":"circulatory system"},
                  {"category":1,"category_string":"congenital anomalies"},
                  {"category":2,"category_string":"dermatologic"},
                  {"category":3,"category_string":"digestive"},
                  {"category":4,"category_string":"endocrine\/metabolic"},
                  {"category":5,"category_string":"genitourinary"},
                  {"category":6,"category_string":"hematopoietic"},
                  {"category":7,"category_string":"infectious diseases"},
                  {"category":8,"category_string":"injuries & poisonings"},
                  {"category":9,"category_string":"mental disorders"},
                  {"category":10,"category_string":"musculoskeletal"},
                  {"category":11,"category_string":"neoplasms"},
                  {"category":12,"category_string":"neurological"},
                  {"category":13,"category_string":"other"},
                  {"category":14,"category_string":"pregnancy complications"},
                  {"category":15,"category_string":"respiratory"},
                  {"category":16,"category_string":"sense organs"},
                  {"category":17,"category_string":"symptoms"}]
)}

async function* _hist_data(hist_args,$,Promises)
{
  let data_for_server = hist_args[0];

    let next_data = await $.ajax({
      url: 'http://127.0.0.1:5000/grab_data',
      dataType: 'json',
      data: JSON.stringify(data_for_server),
      contentType: 'application/json;charset=UTF-8',
      type: 'POST'
    });

    yield Promises.delay(1000, JSON.parse(next_data));
}


async function* _jhist_data(var1_select,var2_select,target_var,$,Promises)
{
  let data_for_server = { cmd:"compute_hist2D", var1: var1_select, var2: var2_select, response:target_var};

    let next_data = await $.ajax({
      url: 'http://127.0.0.1:5000/grab_data',
      dataType: 'json',
      data: JSON.stringify(data_for_server),
      contentType: 'application/json;charset=UTF-8',
      type: 'POST'
    });

    yield Promises.delay(1000, JSON.parse(next_data));
}


async function* _var_comp_stats(var1_select,var2_select,target_var,$,Promises)
{
  let data_for_server = { cmd:"independence_tests", var1: var1_select, var2: var2_select, response:target_var};

    let next_data = await $.ajax({
      url: 'http://127.0.0.1:5000/grab_data',
      dataType: 'json',
      data: JSON.stringify(data_for_server),
      contentType: 'application/json;charset=UTF-8',
      type: 'POST'
    });

    yield Promises.delay(1000, JSON.parse(next_data));
}


function _table_cols(reg_panel_size,point_col)
{
  let base_data = [{name: 'PheCode', width:reg_panel_size*0.07},
                   {name: 'Phenotype', width:reg_panel_size*0.4},
                   {name: 'Count', width:reg_panel_size*0.07},
                   {name:  point_col, width:reg_panel_size*0.07},
                   {name: 'P-value', width:reg_panel_size*0.09},
                   {name: 'Category', width:reg_panel_size*0.15}
                  ]
  let offset = 0;
  for (var i=0; i<base_data.length;i++){
    base_data[i].offset = offset; 
    offset = offset + base_data[i].width
  }
  return base_data
}


function _table_data(z,reg_data,start_ix,num_rows)
{
  let table_data = z.sortByCol('pval','asc',reg_data)
      
  return table_data.slice(start_ix,start_ix+num_rows,table_data)
}


function _86(md){return(
md `### Imports`
)}

function _d3(require){return(
require('d3@5')
)}

function _z(require){return(
require('zebras')
)}

function _math(require){return(
require('https://unpkg.com/mathjs@5.9.0/dist/math.min.js')
)}

export default function define(runtime, observer) {
  const main = runtime.module();
  main.variable(observer()).define(["md"], _1);
  main.variable(observer()).define(["group_panel_size","margin","d3","group_panel"], _2);
  main.variable(observer()).define(["reg_panel_size","margin","d3","reg_panel"], _3);
  main.variable(observer()).define(["md"], _4);
  main.variable(observer("group_panel")).define("group_panel", ["group_panel_size","margin","var_dist_size","legend_height","plot_groupvar_hists","indep_comp","reg_builder","group_legend","d3","var_dist2_size"], _group_panel);
  main.variable(observer("plot_groupvar_hists")).define("plot_groupvar_hists", ["d3","group_vars","var_dist_size","z","hist_data","group_data","target_var","mutable target_var","mutable group_data_args","mutable hist_args","cov_select","mutable cov_select","var_comp","mutable var_comp","var_stat_scale","draw_groupvar_buttons"], _plot_groupvar_hists);
  main.variable(observer("draw_groupvar_buttons")).define("draw_groupvar_buttons", ["var_comp","d3","target_var","z","mutable var_comp","cov_select","mutable cov_select"], _draw_groupvar_buttons);
  main.variable(observer("indep_comp")).define("indep_comp", ["jhist_data","var_dist2_size","z","var1_select","hist_data","var2_select","d3","comp_stats"], _indep_comp);
  main.variable(observer("comp_stats")).define("comp_stats", ["d3","var_comp_stats","var_dist2_size","var1_select","var2_select","z","var_stat_scale","target_var"], _comp_stats);
  main.variable(observer("reg_builder")).define("reg_builder", ["group_data","reg_builder_size","var_dist_size","d3","math","reg_types","mutable reg_types","z","cov_selection"], _reg_builder);
  main.variable(observer("cov_selection")).define("cov_selection", ["z","cov_select","reg_builder_size","d3","gen_display_string","run_status","reg_types","gen_model_string","target_var","reg_list","save_covariate_data","mutable reg_list","mutable reg_args","mutable start_ix","mutable run_status","mutable save_covariate_data","canonical","mutable canonical","hist_data"], _cov_selection);
  main.variable(observer("group_legend")).define("group_legend", ["var_dist2_size","legend_height","d3","target_var","var_stat_scale"], _group_legend);
  main.variable(observer("gen_display_string")).define("gen_display_string", ["canonical","target_var"], _gen_display_string);
  main.variable(observer("gen_model_string")).define("gen_model_string", ["canonical"], _gen_model_string);
  main.variable(observer()).define(["md"], _15);
  main.variable(observer("reg_panel")).define("reg_panel", ["draw_cat_legend","draw_volcano_legend","thresh_select_panel","plot_logOdds","plot_volcano","drawTable","mutable run_status"], _reg_panel);
  main.variable(observer("plot_logOdds")).define("plot_logOdds", ["z","thresh_value","reg_data","d3","reg_plot_size","low_lim_col","up_lim_col","reg_plot_legend_height","buffer","phe_cat_colors","vol_color","xaxis_name","phe_cat_markers","point_col","model_type","mutable model_type"], _plot_logOdds);
  main.variable(observer("plot_volcano")).define("plot_volcano", ["d3","reg_data","reg_plot_size","point_col","reg_plot_legend_height","buffer","phe_cat_colors","vol_color","xaxis_name"], _plot_volcano);
  main.variable(observer("draw_volcano_legend")).define("draw_volcano_legend", ["reg_plot_legend_height","cat_legend_width","d3"], _draw_volcano_legend);
  main.variable(observer("draw_cat_legend")).define("draw_cat_legend", ["cat_legend_width","reg_plot_legend_height","d3","phe_categories","phe_cat_colors","phe_cat_markers"], _draw_cat_legend);
  main.variable(observer("thresh_select_panel")).define("thresh_select_panel", ["reg_panel_size","cat_legend_width","reg_plot_legend_height","d3","thresh_types","mutable thresh_types","z"], _thresh_select_panel);
  main.variable(observer("drawTable")).define("drawTable", ["reg_plot_legend_height","reg_plot_size","reg_panel_size","reg_table_height","d3","reg_data","model_string","table_data","table_cols","z","point_col","phe_cat_colors","vol_color","buffer","start_ix","num_rows","mutable start_ix"], _drawTable);
  main.variable(observer()).define(["md"], _23);
  main.variable(observer("fdr_thresh")).define("fdr_thresh", ["z","reg_data","alpha"], _fdr_thresh);
  main.variable(observer("bon_thresh")).define("bon_thresh", ["z","reg_data","alpha"], _bon_thresh);
  main.variable(observer("thresh_value")).define("thresh_value", ["z","thresh_types","fdr_thresh","bon_thresh"], _thresh_value);
  main.variable(observer()).define(["md"], _27);
  main.variable(observer("phe_cat_colors")).define("phe_cat_colors", _phe_cat_colors);
  main.variable(observer("phe_cat_markers")).define("phe_cat_markers", ["d3"], _phe_cat_markers);
  main.variable(observer("var_stat_scale")).define("var_stat_scale", ["d3"], _var_stat_scale);
  main.variable(observer("vol_color")).define("vol_color", ["bon_thresh","fdr_thresh"], _vol_color);
  main.variable(observer()).define(["md"], _32);
  main.variable(observer("buffer")).define("buffer", _buffer);
  main.variable(observer("margin")).define("margin", _margin);
  main.variable(observer("alpha")).define("alpha", _alpha);
  main.variable(observer()).define(["md"], _36);
  main.variable(observer("reg_panel_size")).define("reg_panel_size", _reg_panel_size);
  main.variable(observer("reg_plot_size")).define("reg_plot_size", ["reg_panel_size","reg_plot_legend_height","reg_table_height"], _reg_plot_size);
  main.variable(observer("reg_plot_legend_height")).define("reg_plot_legend_height", _reg_plot_legend_height);
  main.variable(observer("reg_table_height")).define("reg_table_height", _reg_table_height);
  main.variable(observer("cat_legend_width")).define("cat_legend_width", ["reg_panel_size"], _cat_legend_width);
  main.variable(observer("num_rows")).define("num_rows", _num_rows);
  main.variable(observer("point_col")).define("point_col", ["OR_type","model_type"], _point_col);
  main.variable(observer("up_lim_col")).define("up_lim_col", ["OR_type","model_type"], _up_lim_col);
  main.variable(observer("low_lim_col")).define("low_lim_col", ["z","OR_type","model_type"], _low_lim_col);
  main.variable(observer("xaxis_name")).define("xaxis_name", ["OR_type","model_type"], _xaxis_name);
  main.variable(observer()).define(["md"], _47);
  main.variable(observer("legend_height")).define("legend_height", _legend_height);
  main.variable(observer("group_panel_size")).define("group_panel_size", _group_panel_size);
  main.variable(observer("var_dist_size")).define("var_dist_size", ["group_panel_size"], _var_dist_size);
  main.variable(observer("var_dist2_size")).define("var_dist2_size", ["group_panel_size","legend_height"], _var_dist2_size);
  main.variable(observer("reg_builder_size")).define("reg_builder_size", ["group_panel_size"], _reg_builder_size);
  main.variable(observer("var1_select")).define("var1_select", ["z","var_comp"], _var1_select);
  main.variable(observer("var2_select")).define("var2_select", ["z","var_comp"], _var2_select);
  main.variable(observer()).define(["md"], _55);
  main.variable(observer()).define(["md"], _56);
  main.define("initial target_var", ["response_var"], _target_var);
  main.variable(observer("mutable target_var")).define("mutable target_var", ["Mutable", "initial target_var"], (M, _) => new M(_));
  main.variable(observer("target_var")).define("target_var", ["mutable target_var"], _ => _.generator);
  main.variable(observer("group_vars")).define("group_vars", ["z","group_var_init"], _group_vars);
  main.define("initial group_data_args", _group_data_args);
  main.variable(observer("mutable group_data_args")).define("mutable group_data_args", ["Mutable", "initial group_data_args"], (M, _) => new M(_));
  main.variable(observer("group_data_args")).define("group_data_args", ["mutable group_data_args"], _ => _.generator);
  main.define("initial hist_args", _hist_args);
  main.variable(observer("mutable hist_args")).define("mutable hist_args", ["Mutable", "initial hist_args"], (M, _) => new M(_));
  main.variable(observer("hist_args")).define("hist_args", ["mutable hist_args"], _ => _.generator);
  main.define("initial reg_args", _reg_args);
  main.variable(observer("mutable reg_args")).define("mutable reg_args", ["Mutable", "initial reg_args"], (M, _) => new M(_));
  main.variable(observer("reg_args")).define("reg_args", ["mutable reg_args"], _ => _.generator);
  main.define("initial reg_list", _reg_list);
  main.variable(observer("mutable reg_list")).define("mutable reg_list", ["Mutable", "initial reg_list"], (M, _) => new M(_));
  main.variable(observer("reg_list")).define("reg_list", ["mutable reg_list"], _ => _.generator);
  main.define("initial var_comp", ["group_vars"], _var_comp);
  main.variable(observer("mutable var_comp")).define("mutable var_comp", ["Mutable", "initial var_comp"], (M, _) => new M(_));
  main.variable(observer("var_comp")).define("var_comp", ["mutable var_comp"], _ => _.generator);
  main.define("initial cov_select", ["group_vars"], _cov_select);
  main.variable(observer("mutable cov_select")).define("mutable cov_select", ["Mutable", "initial cov_select"], (M, _) => new M(_));
  main.variable(observer("cov_select")).define("cov_select", ["mutable cov_select"], _ => _.generator);
  main.define("initial reg_types", _reg_types);
  main.variable(observer("mutable reg_types")).define("mutable reg_types", ["Mutable", "initial reg_types"], (M, _) => new M(_));
  main.variable(observer("reg_types")).define("reg_types", ["mutable reg_types"], _ => _.generator);
  main.define("initial thresh_types", _thresh_types);
  main.variable(observer("mutable thresh_types")).define("mutable thresh_types", ["Mutable", "initial thresh_types"], (M, _) => new M(_));
  main.variable(observer("thresh_types")).define("thresh_types", ["mutable thresh_types"], _ => _.generator);
  main.define("initial OR_type", _OR_type);
  main.variable(observer("mutable OR_type")).define("mutable OR_type", ["Mutable", "initial OR_type"], (M, _) => new M(_));
  main.variable(observer("OR_type")).define("OR_type", ["mutable OR_type"], _ => _.generator);
  main.define("initial start_ix", _start_ix);
  main.variable(observer("mutable start_ix")).define("mutable start_ix", ["Mutable", "initial start_ix"], (M, _) => new M(_));
  main.variable(observer("start_ix")).define("start_ix", ["mutable start_ix"], _ => _.generator);
  main.define("initial run_status", _run_status);
  main.variable(observer("mutable run_status")).define("mutable run_status", ["Mutable", "initial run_status"], (M, _) => new M(_));
  main.variable(observer("run_status")).define("run_status", ["mutable run_status"], _ => _.generator);
  main.define("initial save_covariate_data", _save_covariate_data);
  main.variable(observer("mutable save_covariate_data")).define("mutable save_covariate_data", ["Mutable", "initial save_covariate_data"], (M, _) => new M(_));
  main.variable(observer("save_covariate_data")).define("save_covariate_data", ["mutable save_covariate_data"], _ => _.generator);
  main.define("initial canonical", _canonical);
  main.variable(observer("mutable canonical")).define("mutable canonical", ["Mutable", "initial canonical"], (M, _) => new M(_));
  main.variable(observer("canonical")).define("canonical", ["mutable canonical"], _ => _.generator);
  main.define("initial model_type", _model_type);
  main.variable(observer("mutable model_type")).define("mutable model_type", ["Mutable", "initial model_type"], (M, _) => new M(_));
  main.variable(observer("model_type")).define("model_type", ["mutable model_type"], _ => _.generator);
  main.variable(observer("model_string")).define("model_string", ["reg","mutable model_type"], _model_string);
  main.variable(observer()).define(["md"], _74);
  main.variable(observer("response_var")).define("response_var", ["$","Promises"], _response_var);
  main.variable(observer("group_var_init")).define("group_var_init", ["$","Promises"], _group_var_init);
  main.variable(observer("group_data")).define("group_data", ["group_data_args","$","Promises"], _group_data);
  main.variable(observer("reg")).define("reg", ["reg_args","$","Promises"], _reg);
  main.variable(observer("reg_data")).define("reg_data", ["z","reg"], _reg_data);
  main.variable(observer("phe_categories")).define("phe_categories", _phe_categories);
  main.variable(observer("hist_data")).define("hist_data", ["hist_args","$","Promises"], _hist_data);
  main.variable(observer("jhist_data")).define("jhist_data", ["var1_select","var2_select","target_var","$","Promises"], _jhist_data);
  main.variable(observer("var_comp_stats")).define("var_comp_stats", ["var1_select","var2_select","target_var","$","Promises"], _var_comp_stats);
  main.variable(observer("table_cols")).define("table_cols", ["reg_panel_size","point_col"], _table_cols);
  main.variable(observer("table_data")).define("table_data", ["z","reg_data","start_ix","num_rows"], _table_data);
  main.variable(observer()).define(["md"], _86);
  main.variable(observer("d3")).define("d3", ["require"], _d3);
  main.variable(observer("z")).define("z", ["require"], _z);
  main.variable(observer("math")).define("math", ["require"], _math);
  const child1 = runtime.module(define1);
  main.import("jQuery", "$", child1);
  return main;
}
