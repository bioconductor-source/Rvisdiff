var svgStyles = `
text { font-size: 10px; font-family: sans-serif; }

rect.button { fill: #888; stroke: #666; stroke-width: 1.5px; cursor: pointer; }
image.changeColor { cursor: pointer; }
.rowClust .link, .colClust .link { vector-effect: non-scaling-stroke; fill: none; stroke: #ccc; stroke-width: 1.5px; shape-rendering: crispEdges; }
.brush .extent{ fill-opacity: 0.15; stroke-width: 1px; fill: #333; stroke: #333; }

.marker { fill: none; stroke: black; stroke-width: 3px; opacity: 0; }
.axis path, .axis line { fill: none; stroke: #aaa; shape-rendering: crispEdges; }
.axis.grid .tick line { opacity: 0.2; }

path.line { fill: none; }
line.bound { stroke: black; fill: none; stroke-dasharray: 3,3; }

.boxplot line, .boxplot rect, .boxplot circle { stroke: #000; stroke-width: 1px; }
.boxplot rect { fill: #1f77b4; }
.boxplot.g2 rect { fill: #ff7f0e; }
.boxplot .center { stroke-dasharray: 3,3; }
.boxplot .outlier { fill: none;  stroke: #000; }
g.left-arrow, g.right-arrow { fill: #fff;  stroke: #000; stroke-width: 2px; cursor: pointer; }
#boxplot g[clip-path] { cursor: grab; }

.g-a { fill: #aec7e8; }
.g-b { fill: #ffbb78; }
.g-split { stroke: #000; stroke-opacity: .18; shape-rendering: crispEdges; }

.buttons rect { fill: #888; stroke: #666; stroke-width: 1.5px; cursor: pointer; }
.buttons text { stroke-width: 0.5px; font-size: 10px; fill: #444; }
.buttons .axis path, .buttons .axis line { shape-rendering: auto; }
.scale text { font-size: 12px; fill: #444; }
line.link { stroke: #999; stroke-opacity: .6; }
.linkText { font-size: 7px; fill: #999; }
.net > .node > circle { stroke: #fff; stroke-width: 1.5px; }
.node > text { fill: #333; }
.area { opacity: 0.2; stroke-width: 3; }
`;

var groups = false,
    samples = false,
    groupnames = false,
    cutoff = 0.05,
    normalized = false,
    results = false,
    cpms = false,
    rows = false,
    cols = false,
    plotpos = false,
    genes = false,
    idx = {
      genes: 0,
      expMean: null,
      log2FC: null,
      pvalue: null,
      padj: null
    };

window.onload = function(){
  d3.select("head").append("style").text(svgStyles);

  try{
    var json = JSON.parse(d3.select("#data").text());
  }catch(e){
    $(".loading").text("Error");
    console.log(e);
    return;
  }
  cutoff = json.cutoff;
  if(json.normalized){
    normalized = true;
  }
  if(json.groups && json.groups.length){
    groups = json.groups;
    groupnames = groups.filter(function(value, index, self){ return self.indexOf(value) === index; });
  }
  if(json.cpms && json.cpms.length){
    cpms = json.cpms;
    samples = cpms.shift();
  }
  results = json.DE;
  cols = results.shift();
  rows = results.map(function(d){ return d[0]; });
  d3.keys(idx).forEach(function(n,i){
    if(i)
      idx[n] = cols.indexOf(json.names[i-1]);
  });

  d3.select("#information").text(results.filter(function(d){ return d[idx.padj]<cutoff; }).length + " significant genes of " + results.length);

  plotpos = results.map(function(d){ return { gene:d[idx.genes]}; });
  render_plots();
  displayTable($('#table'),results,cols);
  $(".loading").remove();
  $(".content").removeClass('hidden');

  filter_cpms();
  if(cpms)
    $("div.cpmgraphs").removeClass("hidden");
}

window.onresize = function(){
  render_plots();
  filter_cpms(genes);
}

function render_plots(){
  var width = $(window).width()-60,
      height = 400;
  if(width>768){
    width = Math.floor(width/2);
  }

  plot("volcano",results.map(function(d){ return [d[idx.log2FC],-(Math.log(d[idx.pvalue])/Math.log(10))]; }),cols[idx.log2FC],"-log10 p-value",width,height);
  var fnMap = function(d){ return [d[idx.expMean],d[idx.log2FC]]; };
  if(cols[idx.expMean]=="baseMean"){
    fnMap = function(d){ return [Math.log(d[idx.expMean])/Math.log(2),d[idx.log2FC]]; };
  }
  plot("maplot",results.map(fnMap),cols[idx.expMean],cols[idx.log2FC],width,height);
}

function filter_cpms(filter) {
  var highlight = [];

  $("p.caption > .first-genes").addClass("hidden");

  if(!filter || !filter.length) {
    filter = results.filter(function(d){ return d[idx.padj]<cutoff; });
    if(!filter.length)
      filter = results;
    else
      highlight = filter.map(function(d){ return d[idx.genes]; });
  }else
    highlight = filter.map(function(d){ return d[idx.genes]; });

  filter = filter.sort(function(a,b){ return a[idx.pvalue]-b[idx.pvalue]; });
  if(filter.length>100){
    $("p.caption > .first-genes").removeClass("hidden");
    filter = filter.slice(0,100);
  }

  plotpos.forEach(function(d){
    if(highlight.indexOf(d.gene)!=-1)
      d.selected = true;
    else
      delete d.selected;
  })
  plot("volcano");
  plot("maplot");

  if(cpms){

    var indices = filter.map(function(d){ return rows.indexOf(d[idx.genes]); }).filter(function(d){ return d!=-1; });

    var frows = filter.map(function(d){ return d[idx.genes]; }),
        data = indices.map(function(d){ return cpms[d]; }),
        dicGroups = groups ? groups.map(function(d){ return groupnames.indexOf(d); }) : false;
    linechart(frows,samples,data,dicGroups);
    var boxdata = data;
    if(!normalized){
      boxdata = data.map(function(d){ return d.map(function(dd){ return Math.log(dd+1)/Math.log(2); }) });
    }
    boxplot(frows,samples,boxdata,dicGroups);

    filter = filter.sort(function(a,b){ return a[idx.log2FC]-b[idx.log2FC]; });
    indices = filter.map(function(d){ return rows.indexOf(d[idx.genes]); }).filter(function(d){ return d!=-1; });

    frows = filter.map(function(d){ return d[idx.genes]; });
    data = indices.map(function(d){ return cpms[d]; });
    heatmap(frows,samples,data,groups,true);
  }
}

function displayTable(sel, data, cols){

      selectableTable(sel, data, cols, idx.pvalue+1,
        function(table){
          genes = table.rows( { selected: true } ).data();
          filter_cpms(genes);
        });

      sel.find('tbody').on("mouseover", "tr", function(){
        var gene = $(this).find('td:nth-child(2)').text(),
            pos = plotpos.filter(function(d){ return d.gene==gene; })[0];
        ['volcano','maplot'].forEach(function(d){
          var marker = d3.select("div#"+d+" .marker");
          marker
            .attr("cx",pos[d].x)
            .attr("cy",pos[d].y)
            .style("opacity",1)
        })
      }).on("mouseout", "tr", function(){
        d3.selectAll("div.plot .marker")
          .style("opacity",0);
      })

  function selectableTable(sel, data, cols, orderIndex, selectAction, addColumnDef){
      var preparedCols = cols.map(function(d,i){ return {title: d, data: i}; });
      preparedCols.unshift({ title: "", data: null, defaultContent: "" });

      var columnDefs = [{
            orderable: false,
            className: 'select-checkbox',
            targets:   0
        }];
      if(typeof addColumnDef != 'undefined'){
        if(Array.isArray(addColumnDef))
          addColumnDef.forEach(function(d){ columnDefs.push(d); });
        else
          columnDefs.push(addColumnDef);
      }
      columnDefs.push({
            render: function (data, type, row) {
              if (type === 'display') {
                return formatter( data );
              }
              return data;
            },
            targets: '_all'
        });

      var table = sel.DataTable({
        data: data,
        columns: preparedCols,
        columnDefs: columnDefs,
        dom: 'Blftp',
        buttons: [
            {
                text: 'Select all',
                action: function () {
                    table.rows( { page: 'current' } ).select();
                }
            },
            {
                text: 'Select none',
                action: function () {
                    table.rows().deselect();
                }
            }
        ],
        order: [[ orderIndex, "asc" ]],
        select: {
          style: 'multi'
        }
      });

      new $.fn.dataTable.Buttons( table, { buttons: ['copy','csv','excel'] } );
 
      table.buttons( 1, null ).container().appendTo( table.table().container() );

      table.on( 'select', function ( e, dt, type ) {
        trSelection(type);
      }).on( 'deselect', function ( e, dt, type ) {
        trSelection(type);
      });

    function trSelection(type,indexes) {
        if ( type === 'row' ) {
          selectAction(table);
        }
    }
  }
}

function plot(id,data,xLab,yLab,width,height){
  var canvasScale = 4;

  if(typeof data == 'undefined'){
    var selected = [],
        canvas = d3.select("#"+id+" canvas"),
        width = +canvas.attr('width'),
        height = +canvas.attr('height'),
        ctx = canvas.node().getContext("2d");
    ctx.clearRect(0, 0, width, height);

    ctx.fillStyle = "#3182bd";
    ctx.globalAlpha = 0.6;
    plotpos.forEach(function(d, i) {
      if(d.selected)
        selected.push(i);
      else
        drawCircle(ctx,d[id].x,d[id].y);
    })

    ctx.fillStyle = "#B22222";
    ctx.globalAlpha = 0.8;
    selected.forEach(function(i){
      var d = plotpos[i][id];
      drawCircle(ctx,d.x,d.y);
    })
    return;
  }

var margin = {top: 40, right: 40, bottom: 80, left: 90};

  var div = d3.select("#"+id);
  div.selectAll("*").remove();
  div.style("position","relative")
     .style("width",width+"px")
     .style("height",height+"px")

width = width - margin.left - margin.right;
height = height - margin.top - margin.bottom;

var x = d3.scale.linear()
    .range([0, width])
    .domain(d3.extent(data.filter(function(d){ return isFinite(d[0]); }), function(d) { return d[0]; }))
    .nice();

var y = d3.scale.linear()
    .range([height, 0])
    .domain(d3.extent(data.filter(function(d){ return isFinite(d[1]); }), function(d) { return d[1]; }))
    .nice();

data.forEach(function(d){
  if(!isFinite(d[0])){
    if(d[0]>0)
      d[0] = x.domain()[1];
    else
      d[0] = x.domain()[0];
  }
  if(!isFinite(d[1])){
    if(d[1]>0)
      d[1] = y.domain()[1];
    else
      d[1] = y.domain()[0];
  }
})

var xAxis = d3.svg.axis()
    .scale(x)
    .tickFormat(formatter)
    .orient("bottom");

var yAxis = d3.svg.axis()
    .scale(y)
    .tickFormat(formatter)
    .orient("left");

  var tooltip = div.append("div").attr("class","tooltip")

  var svg = div.append("svg")
      .attr("width", width + margin.left + margin.right)
      .attr("height", height + margin.top + margin.bottom)
    .on("mousemove", function(){
      var pos = {
        x: d3.mouse(this)[0],
        y: d3.mouse(this)[1]
      };
      pos.minX = - margin.left + pos.x - 2;
      pos.maxX = - margin.left + pos.x + 2;
      pos.minY = - margin.top + pos.y - 2;
      pos.maxY = - margin.top + pos.y + 2;

      var matches = plotpos.filter(function(d) {
        if (d[id].x >= pos.minX && d[id].x <= pos.maxX && 
            d[id].y >= pos.minY && d[id].y <= pos.maxY) {
          return d;
        }
      })

      if(matches.length){
        var txt = [];
        matches.forEach(function(d){
          txt.push(d.gene+" ("+formatter(x.invert(d[id].x))+","+formatter(y.invert(d[id].y))+")");
        })
        tooltip.html(txt.join("<br/>"));
        tooltip.style({"display": "block",
          "left": (pos.x + 10) + "px",
          "top" : (pos.y + 10) + "px"});
      }else{
        tooltip.html("");
        tooltip.style("display", "none");
      }
    })

  var canvas = svg.append("foreignObject")
     .attr({"x":margin.left,"y":margin.top,"width":width,"height":height})
     .append("xhtml:body")
       .style({"margin":"0px","padding":"0px","background-color":"none","width":width+"px","height":height+"px"})
       .append("canvas")
         .attr("width", width*canvasScale)
         .attr("height", height*canvasScale)
         .style("width", width+"px")
         .style("height", height+"px");

  var ctx = canvas.node().getContext("2d");

  ctx.clearRect(0, 0, width*canvasScale, height*canvasScale);

  ctx.fillStyle = "#3182bd";
  ctx.globalAlpha = 0.6;

  var dx, dy;
  data.forEach(function(d, i) {
      dx = x(d[0]);
      dy = y(d[1]);
      plotpos[i][id] = {x:dx, y:dy};
      drawCircle(ctx,dx,dy);
  })

  svg = svg.append("g")
      .attr("transform", "translate(" + margin.left + "," + margin.top + ")");

  svg.append("g")
      .attr("class", "x axis")
      .attr("transform", "translate(0," + height + ")")
      .call(xAxis)
      .append("text")
        .attr("x", width/2)
        .attr("y", margin.bottom-10)
        .style("text-anchor", "middle")
        .text(xLab);

  svg.append("g")
      .attr("class", "y axis")
      .call(yAxis)
      .append("text")
        .attr("transform", "rotate(-90)")
        .attr("y", -margin.left+30)
        .attr("x", -height/2)
        .style("text-anchor", "middle")
        .text(yLab)

  svg.append("circle")
      .attr("class", "marker")
      .attr("r",6)

  function drawCircle(ctx,dx,dy){
      ctx.beginPath();
      ctx.arc(dx*canvasScale, dy*canvasScale, 3*canvasScale, 0, 2 * Math.PI, true);
      ctx.closePath();
      ctx.fill();
  }

  pngExportButton(div,['png'],{ volcano: "VolcanoPlot", maplot: "MAplot" }[id]);
}

function linechart(rows,cols,data,groups){
var width = $(window).width() - 60,
    height = 400,
    margin = {top: 40, right: 40, bottom: 80, left: 90},
    xLab = "samples",
    yLab = "expression";

if(width>768)
  width = Math.floor(width/2);

width = width - margin.left - margin.right;
height = height - margin.top - margin.bottom;

var color = d3.scale.category10();

var x = d3.scale.ordinal()
    .rangePoints([0, width])
    .domain(d3.range(0,cols.length));

var y = d3.scale.linear()
    .range([height, 0])
    .domain([0,d3.max([].concat.apply([], data))])
    .nice();

var xAxis = d3.svg.axis()
    .scale(x)
    .tickFormat(function(d){ return cols[d]; })
    .orient("bottom");

var yAxis = d3.svg.axis()
    .scale(y)
    .tickFormat(formatter)
    .orient("left");

var lineFunction = d3.svg.line()
      .x(function(d,i) { return x(i); })
      .y(function(d) { return y(d); })
      .interpolate("linear");

  var div = d3.select("#linechart");

  div.selectAll("*").remove();

  var tooltip = div.append("div").attr("class","tooltip")

  var svg = div.append("svg")
      .attr("width", width + margin.left + margin.right)
      .attr("height", height + margin.top + margin.bottom);

  var defs = svg.append("defs");

  defs.append("clipPath")
      .attr("id","lines-clip")
      .append("path")
        .attr("d","M0,0h"+width+"v"+height+"h"+(-width)+"Z")

  var g = svg.append("g")
      .attr("transform", "translate(" + margin.left + "," + margin.top + ")");

  g.append("g")
      .attr("class", "x axis")
      .attr("transform", "translate(0," + (height+10) + ")")
      .call(xAxis)
      .select(".domain").remove()

  g.selectAll(".x.axis .tick text")
	  .attr("x", 8)
	  .attr("y", 8)
	  .attr("transform", "rotate(45)")
	  .style("text-anchor", "start")

  var gyAxis = g.append("g")
      .attr("class", "y axis")
      .call(yAxis)
  gyAxis.append("text")
        .attr("transform", "rotate(-90)")
        .attr("y", -margin.left+30)
        .attr("x", -height/2)
        .style("text-anchor", "middle")
        .text(yLab)

  if(groups){
    var bound = groups.indexOf(1);
    bound = (x(bound) + x(bound - 1)) / 2;
    g.append("line")
      .attr("class","bound")
      .attr("x1",bound)
      .attr("y1",0)
      .attr("x2",bound)
      .attr("y2",height)
  }

  var lines = g.append("g")
                 .attr("clip-path", "url(#lines-clip)")
        .selectAll(".line")
      .data(data)
    .enter().append("path")
      .attr("class", "line")
      .attr("d",lineFunction)
      .style("stroke", function(d,i){ return color(i); })
      .on("mouseover", function(d,i){
          d3.select(this).style("stroke-width","3px");
          tooltip.text(rows[i]);
          tooltip.style({"display": "block",
            "left": d3.mouse(div.node())[0] + 10 + "px",
            "top" : d3.mouse(div.node())[1] + 10 + "px"});
      })
      .on("mouseout", function(){
          d3.select(this).style("stroke-width",null);
          tooltip.text("");
          tooltip.style("display","none");
      });

  brushSlider(svg,y.domain(),function(val){
    y.domain(val)
    gyAxis.call(yAxis)
    lines.attr("d",lineFunction)
  })

  pngExportButton(div,['png','svg'],"LinesDiagram");
}

function boxplot(rows,cols,data,groups){
var width = $(window).width() - 60,
    height = 400,
    margin = {top: 40, right: 40, bottom: 80, left: 90},
    xLab = "genes",
    yLab = "expression",
    maxData = d3.max([].concat.apply([], data)),
    minData = d3.min([].concat.apply([], data));

if(width>768)
  width = Math.floor(width/2);

width = width - margin.left - margin.right;
height = height - margin.top - margin.bottom;

var color = d3.scale.category10();

  var x = d3.scale.ordinal()
                .domain(rows)
                .rangeBands([0, (width/5)*rows.length ], 0.6, 0.3);

  var xAxis = d3.svg.axis()
                .scale(x)
                .tickFormat(function(d){ return d; })
                .orient("bottom");

  var y = d3.scale.linear()
                .range([height, 0])
                .domain([minData,maxData])
                .nice();

  var yAxis = d3.svg.axis()
                .scale(y)
                .tickFormat(formatter)
                .orient("left");

  var chart = d3.box()
                .whiskers(iqr(1.5))
                .height(height)
                .width(x.rangeBand()/2)
                .domain(y.domain())
                .tickFormat(d3.format(".1f"))

  var div = d3.select("#boxplot");

  div.selectAll("*").remove();

  var svg = div.append("svg")
      .attr("width", width + margin.left + margin.right)
      .attr("height", height + margin.top + margin.bottom);

  var defs = svg.append("defs");

  defs.append("clipPath")
      .attr("id","boxplot-clip")
      .append("path")
        .attr("d","M"+margin.left+",0h"+width+"v"+(margin.top+height+margin.bottom)+"h"+(-width)+"Z")

  defs.append("clipPath")
      .attr("id","gene-clip")
      .append("path")
        .attr("d","M"+(-x.rangeBand())+",-4h"+x.rangeBand()*3+"v"+(height+8)+"h"+(-x.rangeBand()*3)+"Z")

  var gyAxis = svg.append("g")
      .attr("class", "y axis")
      .attr("transform", "translate(" + margin.left + "," + margin.top + ")")
      .call(yAxis)

  gyAxis.append("text")
        .attr("transform", "rotate(-90)")
        .attr("y", -margin.left+30)
        .attr("x", -height/2)
        .style("text-anchor", "middle")
        .text(yLab)

  var sliderparent = svg.append("g")
        .attr("clip-path", "url(#boxplot-clip)")

  sliderparent.append("rect")
  .attr("fill","transparent")
  .attr("width",width)
  .attr("height",height)
  .attr("x",margin.left)
  .attr("y",margin.top)

  var slider = sliderparent.append("g")
          .attr("transform", "translate(" + margin.left + "," + margin.top + ")");

  var boxes = slider.selectAll("g.gene")	   
        .data(rows)
      .enter().append("g")
        .attr("class","gene")
        .attr("clip-path", "url(#gene-clip)")
        .attr("transform", function(d) { return "translate(" +  x(d)  + ",0)"; } )
        .each(function(d,i){
    if(groups){
      d3.select(this).append("g")
        .datum(data[i].filter(function(d,i){ return groups[i]==0; }))
        .attr("class","boxplot g1")
        .attr("transform", function(d,i) { return "translate(" +  -(x.rangeBand()/4)  + ",0)"; } )
        .call(chart);
      d3.select(this).append("g")
        .datum(data[i].filter(function(d,i){ return groups[i]==1; }))
        .attr("class","boxplot g2")
        .attr("transform", function(d,i) { return "translate(" +  (x.rangeBand()/4)*3  + ",0)"; } )
        .call(chart);
    }else{
      d3.select(this).append("g")
        .datum(data[i])
        .attr("class","boxplot g1")
        .attr("transform", function(d,i) { return "translate(" +  (x.rangeBand()/4)  + ",0)"; } )
        .call(chart);
    }

        })

  boxes.selectAll("text").style("opacity",0);
  boxes.selectAll(".boxplot")
    .on("mouseover",function(){
      d3.select(this).selectAll("text")
          .transition()
          .duration(250)
          .style("opacity",1);
    })
    .on("mouseout",function(){
      d3.select(this).selectAll("text")
          .transition()
          .duration(250)
          .style("opacity",0);
    })

  slider.append("g")
      .attr("class", "x axis")
      .attr("transform", "translate(0," + (height+10) + ")")
      .call(xAxis)
      .select(".domain").remove()

  slider.selectAll(".x.axis .tick text")
	  .attr("x", 8)
	  .attr("y", 8)
	  .attr("transform", "rotate(45)")
	  .style("text-anchor", "start")

  brushSlider(svg,[minData,maxData],function(val){
    y.domain(val)
    chart.domain(val)
    gyAxis.call(yAxis)
    boxes.selectAll("g.boxplot").call(chart)
    boxes.selectAll("text").style("opacity",0);
  });
  arrow(svg,0);
  arrow(svg,1);

  var pos;
  sliderparent
  .on("mousedown",function(){
    sliderparent.style("cursor","grabbing");
    pos = d3.mouse(this)[0];
  })
  .on("mouseup",function(){
    sliderparent.style("cursor",null);
    pos = pos - d3.mouse(this)[0];
    if(pos!=0){
      if(pos > 0){
        pos = 1;
      }else{
        pos = 0;
      }
      moveBoxes(pos);
    }
  })

  function arrow(svg,dir){
    var g = svg.append("g")
       .attr("class",(dir?"right":"left")+"-arrow")
       .attr("transform","translate(" + (dir? margin.left+width+(margin.right/2) : margin.left/2) + "," + (margin.top+height-6) + ")")
       .on("click",function(){
         moveBoxes(dir);
       })

       g.append("circle")
       .attr("cx",0)
       .attr("cy",0)
       .attr("r",6)

       g.append("path")
       .attr("d",dir?"M-2,-3L2,0L-2,3":"M2,-3L-2,0L2,3")
  }

  function moveBoxes(dir){
    var pos = d3.transform(slider.attr("transform")).translate,
        newx = dir? pos[0]-width : pos[0]+width;
    if((newx-margin.left)%width==0 && newx<=margin.left && newx>margin.left-((width*Math.ceil(rows.length/5))))
      slider.transition()
            .duration(750)
            .attr("transform","translate("+newx+","+pos[1]+")")
  }

  function iqr(k) {
    return function(d, i) {
      var q1 = d.quartiles[0],
          q3 = d.quartiles[2],
          iqr = (q3 - q1) * k,
          i = -1,
          j = d.length;
      while (d[++i] < q1 - iqr);
      while (d[--j] > q3 + iqr);
      return [i, j];
    };
  }

  pngExportButton(div,['png','svg'],"Boxplots");
}


function brushSlider(svg,prevVal,callback){

    var y = d3.scale.linear()
            .domain(prevVal)
            .range([100, 0])
            .clamp(true);

    var brush = d3.svg.brush()
            .y(y)
            .extent(prevVal)
            .on("brush", brushed);

    var slider = svg.append("g")
            .attr("class", "brushSlider")
            .attr("transform", "translate(20, 20)");

    slider.append("g")
            .attr("class", "y axis")
            .call(d3.svg.axis()
                    .scale(y)
                    .orient("right")
                    .tickFormat("")
                    .tickSize(0))
            .select(".domain")
            .select(function() { return this.parentNode.appendChild(this.cloneNode(true)); })
            .attr("class", "halo");

    var ruler = slider.append("g")
            .attr("transform", "translate(-4, 0)")
            .attr("class", "ruler")
            .call(brush);

    ruler.selectAll("rect")
        .attr("width", 10);

    ruler.selectAll(".resize > rect")
      .style("visibility","visible")
      .attr("rx",2)
      .attr("x",-1)

    ruler
      .call(brush.extent(prevVal))
      .call(brush.event)

    function brushed() {
        var values = brush.extent();

        if(values[0]==values[1]){
          values = prevVal;
          ruler.call(brush.extent(prevVal))
        }
        callback(values);
    }
}

function heatmap(rows,cols,data,metadata,scaled){
  var width = 0,
      height = 0,
      margin = {top: 70, right: 200, bottom: 140, left: 110},
      cellSize = 14,
      colorScale = scaled?"RdBkGn":"GrRd",
      divHeatmap = d3.select("#heatmap");
  divHeatmap.style("height",divHeatmap.style("height"));
  divHeatmap.selectAll("*").remove();

  width = $(window).width() - 60 - margin.left - margin.right;

 var NAcolor = "transparent";

  var colorScales = {
        GrRd: ["#eeeeee","#E68E8A","#de2d26"],
        GrGn: ["#eeeeee","#90C9A1","#31a354"],
        GrBl: ["#eeeeee","#90B8D6","#3182bd"],
        RdBkGn: ["#de2d26","#000000","#31a354"],
        RdWhBu: ["#de2d26","#ffffff","#3182bd"]
    },
    colorPNG = 'data:image/png;base64,iVBORw0KGgoAAAANSUhEUgAAABQAAAAUCAIAAAAC64paAAADcklEQVQ4y22UTWhcZRiFz/t+3/2ZuTPJ5JqkmcQ0TWiaosVqunBRxIglCF25FelCcOkiBVeiUMR9t6Ibt4KuREQs6s4gahG6CVXaRFrzN5P7M/fnu9+PiwFJ0bM+z/vCgXMoTVOcErGn6qPB0Z00u5upw5HxSfbj7qWV/stREBvbPGE+BZPRxd7OrST9uDtlyEfqcKAwqMWoag9PotXpG5vPbYVeB3BPwEScD+893L4208/CKYg2asZA45HCUY20QVLgOIuTw/mbG58tTa85ZwHw+EQ+vPfop6vLs9lshK6HQEAKCAFiOIGGMZJI/UHR++uD7954eLQz/swAjC53f9w8G+tJH6GDZ8AG5OAsrEVjUWikGrlDLpNmorz17U2AADCx9/iXD+c7SYfhO7CBUBAKpOEcDFBoJA1ShZFGYaKha+0zf/T1bSZmXQ6q+59MOvgKVIFrcAWuIRrAQGlkDYYVhpVIqv6wWD4s5/ab3p0Hf0hmebL7fS+sfQOuAQNokA8BMOAYpcKgxEExcVT0k2q6aCYL1SXdAvmf//yDrA5/78EJDViAAAmSEDW4hBHIKjrOVo6zuZMqzlWnNm1tfbbCOt6+vyNNti8q8JgcR08ghmCQ9Wz9gi7PmTqGiXzrSwg4Qc61LZqKJRRRDggAgAMsYEENhKJWdWVSXZ6xT+fU8aVXMxsmR0RwbTShiqQM+qaAZfCYVIACFfDylWh0YV4vn4gZGQSHASc+V5IaCUc2oHqCQxmcuVTd9QwaaYAGqIASyFucr4bFYmymFr1Q+DIKxDCkNKDcJyUNs1pfOyd7F1958EUYe41fgyqgBkpgdIFGS1JNt1zQc1DWCeMiTYnGiXaj0Bhfbb70LMv2VPD822kpmwIogAqon0K9ArVApiPAobMTxk5qF9dutnD9wswq9eLqDAB2pulff/fvej5roDWgGGoNaolsTC5gsA8K4drWhc5FxkxYFVP15lsb1joGIILOwtaXe3omoXZjz0OtQPedaRMEAQLOdzaA8aBJKiOy17c2mPnfVrlo8eLce9/sBusH/EyGBYWuJqnJWbKWLIQhqWxYqonR9fdf7Z8/898xAEB/3v6UtvNuFXs2rMEp21zYkdCZ30RX5q++c40F/++SjDtKzOLxV78lv+41lS7J6ZaYXj+7+tplAM66095/AIXCz8hIhvLCAAAAAElFTkSuQmCC';

  if(metadata)
    margin.top = margin.top + cellSize + 4;

  height = rows.length*cellSize;

  var collapsed = [];
  if(scaled){
    var mean, sd, i, j;
    for(i = 0; i<data.length; i++){
      mean = d3.mean(data[i]);
      sd = d3.deviation(data[i]);
      for(j = 0; j<data[i].length; j++){
        collapsed.push((data[i][j]-mean)/sd);
      }
    }
  }else
    collapsed = [].concat.apply([], data);

  d3.select("body")
    .on("click",function(){
      d3.select(this).selectAll(".scalePicker, .export-menu").remove();
    })

  var tooltip = divHeatmap.append("div").attr("class","tooltip")

  var svg = divHeatmap.append("svg")
      .attr("width", width + margin.left + margin.right)
      .attr("height", height + margin.top + margin.bottom);

  var defs = svg.append("defs");

  addGradient(defs,'gradButton',['#080','#8f0','#080']);
  d3.entries(colorScales).forEach(function(d){ addGradient(defs,d.key,d.value); });

  var rowClust = svg.append("g")
      .attr("class","rowClust")
      .attr("clip-path",addMask(defs,"rowClust",height,margin.left))
      .attr("transform", "translate(0,"+margin.top+")rotate(90)scale(1,-1)");
  var colClust = svg.append("g")
      .attr("class","colClust")
      .attr("clip-path",addMask(defs,"colClust",width,margin.top))
      .attr("transform", "translate("+margin.left+",0)");
  var gMatrix = svg.append("g")
      .attr("class","matrix")
      .attr("clip-path",addMask(defs,"matrix",width,height))
      .attr("transform", "translate("+margin.left+","+margin.top+")");
  var rowNames = svg.append("g")
      .attr("class","rowNames axis")
      .attr("transform", "translate("+(width+margin.left+4)+","+margin.top+")");
  var colNames = svg.append("g")
      .attr("class","colNames axis")
      .attr("transform", "translate("+margin.left+","+(height+margin.top+4)+")");
  var buttons = svg.append("g")
      .attr("class","buttons")
      .attr("transform", "translate("+(margin.left+width)+",10)");

  displayButtons(buttons);

  dendrogram(rowClust,false,height,margin.left-4);
  dendrogram(colClust,false,width,(metadata?70:margin.top)-4);

  drawMetadata(colClust,metadata,cols);
  drawMatrix(gMatrix,{data:collapsed, rows:rows, cols:cols},colorScale);

  displayNames(rowNames,rows,[0,height],"right");
  displayNames(colNames,cols,[0,width],"bottom");

  divHeatmap.style("height",null);

  function displayButtons(svg){
    svg.append("rect")
      .attr("class","button")
      .attr("x",14)
      .attr("rx",2)
      .attr("width",30)
      .attr("height",10)
      .style("fill",scaled?"url(#gradButton)":"#888")
      .on("click",function(){
        heatmap(rows,cols,data,metadata,scaled?false:true);
      })

    svg.append("text")
      .attr("x",54)
      .attr("y",8)
      .text("rescale")
  }

  function dendrogram(svg,root,width,height){
    if(root){
      var cluster = d3.layout.cluster()
        .separation(function() { return 1; })
        .size([width, height]);

      var diagonal = (function(d) { return "M"+d.source.x+","+d.source.y+"L"+d.target.x+","+d.source.y+"L"+d.target.x+","+d.target.y; });

      var y = d3.scale.linear()
        .range([0, height])

      var nodes = cluster.nodes(root);
          y.domain(d3.extent(nodes, function(d){ return d.height; }).reverse());
          nodes.forEach(function(d){
	    d.y = d.children? y(parseFloat(d.height)) : height;
	  });
      var links = cluster.links(nodes);

      svg.append("g")
      .selectAll(".link")
        .data(links)
      .enter().append("path")
        .attr("class", "link")
        .attr("d", diagonal)
    }else
      svg.append('g')
  }

  function drawMetadata(svg,metadata,cols){
    if(metadata){

      var rows = ["groups"],
          x = d3.scale.linear().domain([0, cols.length]).range([0, width]),
          y = d3.scale.linear().domain([0, rows.length]).range([70, margin.top-4]);

      var cellWidth = x(1)-x(0),
          cellHeight = y(1)-y(0);

      var color = d3.scale.category10(),
          j = 0;

      svg.select("g").selectAll("."+rows[j])
        .data(metadata)
      .enter().append("rect")
        .attr("class","metacell "+rows[j])
        .attr("width", cellWidth)
        .attr("height", cellHeight)
        .attr("x", function(d,i){ return x(i); })
        .attr("y", function(d,i){ return y(j); })
        .style("fill", function(d,i) { return color(d); })
        .append("title")
          .text(function(d,i){ return (rows[j])+"\nsample: "+(cols[i])+"\nvalue: "+d; });
    }
  }

  function drawMatrix(gMatrix,matrix,color){

    var datamax = d3.max(matrix.data),
        datamin = d3.min(matrix.data),
        colorDomain = [datamin, (datamin+datamax)/2, datamax];

    displayScale(colorDomain);

    var rows = matrix.rows.length,
        cols = matrix.cols.length,
        x = d3.scale.linear().domain([0, cols]).range([0, width]),
        y = d3.scale.linear().domain([0, rows]).range([0, height]),
        c = d3.scale.linear().range(colorScales[color]).domain(colorDomain);

    var cell = gMatrix.selectAll(".cell")
      .data(matrix.data)
    .enter().append("rect")
      .attr("class","cell")
      .call(drawCells);
    fillCells(color);

    var brush = d3.svg.brush()
        .x(x)
        .y(y)
        .clamp([true, true])
        .on('brush', function() {
          var extent = brush.extent();
          extent[0][0] = Math.round(extent[0][0]);
          extent[0][1] = Math.round(extent[0][1]);
          extent[1][0] = Math.round(extent[1][0]);
          extent[1][1] = Math.round(extent[1][1]);
          d3.select(this).call(brush.extent(extent));
        })
        .on('brushend', function() {
          zoom(brush.extent());
          brush.clear();
          d3.select(this).call(brush);
        });

    gMatrix.append("g")
      .attr("class", "brush buttons")
      .call(brush)
      .select("rect.background")
        .on("mouseenter",function(){ 
          tooltip.style("display","block");
        })
        .on("mousemove",function(){          
          var col = Math.floor(x.invert(d3.mouse(this)[0]));
          var row = Math.floor(y.invert(d3.mouse(this)[1]));
          var label = formatter(matrix.data[row*cols + col]);
          tooltip.html("row: "+matrix.rows[row]+"<br/>col: "+matrix.cols[col]+"<br/>value: "+label);
          tooltip.style({"left":d3.mouse(divHeatmap.node())[0]+10+"px","top":d3.mouse(divHeatmap.node())[1]+10+"px"});
        })
        .on("mouseleave",function(){ 
          tooltip.style("display","none");
        })

    iconButton(divHeatmap,"colors",colorPNG,"Change Color Scale",scalePicker,{"position":"absolute","top":"30px","left":buttonLeftPos(divHeatmap)});
    pngExportButton(divHeatmap,['png','svg'],"Heatmap");

    function drawCells(cells){
      var cellWidth = x(1)-x(0),
          cellHeight = y(1)-y(0);
      cells
      .attr("width", cellWidth)
      .attr("height", cellHeight)
      .attr("x", function(d,i){ return x(i%cols); })
      .attr("y", function(d,i){ return y(Math.floor(i/cols)); })
    }

    function fillCells(color){
      c.range(colorScales[color]);
      d3.select(".scale rect").attr("fill", "url(#"+color+")")
      d3.selectAll(".cell").transition().duration(500).style("fill", function(d,i) {
        var val = matrix.data[i];
        if(val == null)
          return NAcolor;
        return c(val);
      })
    }

    function displayScale(domain){
      var scaleWidth = 60;

      var x = d3.scale.linear()
      .range([0,scaleWidth/2])
      .domain(domain);

      var axis = d3.svg.axis()
      .scale(x)
      .tickValues(domain)
      .orient("bottom");

      var scale = svg.append("g")
	.attr("class","scale")
        .attr("transform", "translate(20,10)");
        scale.append("rect")
	.attr({x:0, y:0, height:10, width:scaleWidth, rx:2, fill:"black"});
        scale.append("g")
        .attr("class","axis")
        .attr("transform","translate(0,13)")
        .call(axis)
          .select("path.domain").remove()
    }

    function scalePicker(){
      d3.event.stopPropagation();
      if(!d3.select(".scalePicker").empty())
        return;

      var colors = d3.keys(colorScales);

      var picker = svg.append("g")
      .attr("class","scalePicker")
      .attr("transform","translate("+(margin.left + width + margin.right - 90)+",30)");

      picker.append("rect")
      .attr("x",0)
      .attr("y",0)
      .attr("rx",2)
      .attr("width",60)
      .attr("height", 8 + colors.length*14)
      .style({"fill":"white","stroke":"#ccc"})

      picker.selectAll("rect>rect")
        .data(colors)
      .enter().append("rect")
      .attr("x",10)
      .attr("y",function(d,i){ return 6 + i*14; })
      .attr("rx",2)
      .attr("width",40)
      .attr("height",10)
      .attr("fill",function(d){ return "url(#"+d+")"; })
      .style("cursor","pointer")
      .on("click",fillCells);
    }

    function zoom(ex){
      var scale = [1,1], translate = [0,0];

      if(ex[0][0]==ex[1][0] || ex[0][1]==ex[1][1]){
        height = rows * cellSize;
        ex = [[0,0],[cols,rows]];
      }else{
        height = (ex[1][1] - ex[0][1]) * cellSize;
        scale = [
          cols / (ex[1][0] - ex[0][0]),
          rows / (ex[1][1] - ex[0][1])
        ];
        translate = [
          ex[0][0] * (width / cols) * scale[0] * -1,
          ex[0][1] * (height / rows) * scale[1] * -1
        ];
      }

      x.range([translate[0], width * scale[0] + translate[0]]);
      y.range([translate[1], height * scale[1] + translate[1]]);
      drawCells(cell.transition().duration(500));
      displayNames(d3.select("g.rowNames").transition().duration(500),matrix.rows.slice(ex[0][1],ex[1][1]),[0,height],"right");
      displayNames(d3.select("g.colNames").transition().duration(500).attr("transform","translate("+margin.left+","+(height+margin.top+4)+")"),matrix.cols.slice(ex[0][0],ex[1][0]),[0,width],"bottom");
      d3.select(".rowClust>g").transition().duration(500).attr("transform","translate("+translate[1]+",0)scale(1,1)");
      d3.select(".colClust>g").transition().duration(500).attr("transform","translate("+translate[0]+",0)scale("+scale[0]+",1)");
      d3.select("#rowClustMask>rect").transition().duration(500).attr("width",height);
      d3.select("#matrixMask>rect").transition().duration(500).attr("height",height);
      svg.transition().duration(500).attr("height", height + margin.top + margin.bottom);
    }
  }

  function displayNames(svg,names,range,orient){

    var x = d3.scale.ordinal()
    .rangeBands(range)
    .domain(names);

    var xAxis = d3.svg.axis()
    .scale(x)
    .tickFormat(String)
    .orient(orient);

    svg.call(xAxis).select("path.domain").remove();

    if(orient == "bottom")
      svg.selectAll(".tick text")
	.attr("x", 6)
	.attr("y", 6)
	.attr("transform", "rotate(45)")
	.style("text-anchor","start")
  }

  function addGradient(defs, id, stops){
    var offset = 100/(stops.length-1);
    var gradient = defs.append("linearGradient")
	.attr("id",id)
	.attr("x1","0%")
	.attr("y1","0%")
	.attr("x2","100%")
	.attr("y2","0%");

    stops.forEach(function(d, i){
      gradient
	.append("stop")
	.attr("offset",(offset*i)+"%")
	.style("stop-color",d);
    })
  }

  function addMask(defs,id,w,h){
    defs.append("clipPath")
    .attr("id", id+"Mask")
    .append("rect")
        .attr("x", 0)
        .attr("y", 0)
        .attr("width", w)
        .attr("height", h);

    return "url(#"+id+"Mask)";
  }
}

function pngExportButton(sel,options,title){
  var job = null;
  if(options.length==1){
    if(options[0]=="png"){
      job = function(){
        getPNG(sel.select("svg").node(),title);
      };
    }else if(options[0]=="svg"){
      job = function(){
        getSVG(sel.select("svg").node(),title);
      };
    }
  }else if(options.length>1){
    job = function(){
      d3.event.stopPropagation();
      var body = d3.select("body");
      body.selectAll(".export-menu").remove();
      var coor = d3.mouse(body.node());
      var menu = body.append("ul")
      .attr("class","export-menu")
      .style("top",(coor[1]+10)+"px")
      .style("left",(coor[0]-25)+"px")
      options.forEach(function(opt){
        if(opt=="png"){
          menu.append("li").text("PNG").on("click",function(){
            getPNG(sel.select("svg").node(),title);
          })
        }else if(opt=="svg"){
          menu.append("li").text("SVG").on("click",function(){
            getSVG(sel.select("svg").node(),title);
          })
        }
      })
    };
  }
  var svgIcon = "data:image/svg+xml;base64,PD94bWwgdmVyc2lvbj0iMS4wIiBlbmNvZGluZz0iVVRGLTgiIHN0YW5kYWxvbmU9Im5vIj8+CjxzdmcgeG1sbnM6cmRmPSJodHRwOi8vd3d3LnczLm9yZy8xOTk5LzAyLzIyLXJkZi1zeW50YXgtbnMjIiB4bWxucz0iaHR0cDovL3d3dy53My5vcmcvMjAwMC9zdmciIGhlaWdodD0iMTQiIHdpZHRoPSIxNCIgdmVyc2lvbj0iMS4xIiB4bWxuczpjYz0iaHR0cDovL2NyZWF0aXZlY29tbW9ucy5vcmcvbnMjIiB4bWxuczpkYz0iaHR0cDovL3B1cmwub3JnL2RjL2VsZW1lbnRzLzEuMS8iPgogPGVsbGlwc2Ugcng9IjYuNTI4IiBzdHJva2U9IiNjMGMwYzAiIHJ5PSI2LjUyOCIgY3k9IjciIGN4PSI3IiBzdHJva2Utd2lkdGg9Ii45NDQiIGZpbGw9IiNlMGUwZTAiLz4KIDxwYXRoIGQ9Im00LjI4NTIgMi4zNjY2djQuMzI1NWgtMi43OTgzbDIuNzU2MSAyLjg3MDkgMi43NTcgMi44NzEgMi43NTctMi44NzEgMi43NTYtMi44NzA5aC0yLjc5ODN2LTQuMzI1NWgtNS40MjkzeiIgZmlsbD0iIzYwNjA2MCIvPgo8L3N2Zz4K";
  iconButton(sel,"png",svgIcon,"Export",job,{"position":"absolute","top":"10px","left":buttonLeftPos(sel)})
}

function getSVGstring(svg){
  var sheets = document.styleSheets;
  var styleStr = '.buttons,.brushSlider,.left-arrow,.right-arrow';
  Array.prototype.forEach.call(sheets, function(sheet){
    try{ // we need a try-catch block for external stylesheets that could be there...
        styleStr += Array.prototype.reduce.call(sheet.cssRules, function(a, b){
          return a + b.cssText;
        }, "");
    }catch(e){
      //console.log(e);
    }
  });
  var style = document.createElementNS('http://www.w3.org/2000/svg', 'style');
  style.innerHTML = styleStr;
  svg.insertBefore(style, svg.firstElementChild);

  var svgString = new XMLSerializer().serializeToString(svg);
  svg.removeChild(style);

  return svgString;
}

function getPNG(svg,title){
  if(!title){
    title = "graph";
  }

  var svgWidth = svg.getAttribute("width"),
      svgHeight = svg.getAttribute("height"),
      ratio = ratio = 4000 / svgWidth;

  var canvas = document.createElement("canvas");
  canvas.width = svgWidth*ratio;
  canvas.height = svgHeight*ratio;
  var ctx = canvas.getContext("2d");
  ctx.fillStyle = "#ffffff";
  ctx.fillRect(0, 0, canvas.width, canvas.height);
  var innerCanvas = svg.querySelector("canvas");
  if(innerCanvas){
    var x = innerCanvas.parentNode.parentNode.getAttribute("x"),
        y = innerCanvas.parentNode.parentNode.getAttribute("y"),
        w = parseInt(innerCanvas.style.width),
        h = parseInt(innerCanvas.style.height);
    ctx.drawImage(innerCanvas, x*ratio, y*ratio, w*ratio, h*ratio);
    svg.querySelector('foreignObject').remove();
  }

  var styles = document.createElement("style");
  styles.textContent = svgStyles;
  svg.appendChild(styles);
  svg.querySelectorAll(".brushSlider, .left-arrow, .right-arrow").forEach(function(e){ e.remove(); });
  svg.setAttribute("viewBox","0 0 "+svgWidth+" "+svgHeight);
  svg.setAttribute("width",svgWidth*ratio);
  svg.setAttribute("height",svgHeight*ratio);

  var DOMURL = self.URL || self.webkitURL || self;
  var img = new Image();
  var svgString = getSVGstring(svg);
  var svg = new Blob([svgString], {type: "image/svg+xml;charset=utf-8"});
  var url = DOMURL.createObjectURL(svg);
  img.onload = function() {
    ctx.drawImage(img, 0, 0);
    canvas.toBlob(function(blob){
      fileDownload(blob, title+'.png');
    })
  };
  img.src = url;
  window.dispatchEvent(new Event('resize'));
}

function getSVG(svg,title){
  if(!title){
    title = "graph";
  }

  var styles = document.createElement("style");
  styles.textContent = svgStyles;
  svg.appendChild(styles);
  svg.querySelectorAll(".brushSlider, .left-arrow, .right-arrow").forEach(function(e){ e.remove(); });

  var svgString = new XMLSerializer().serializeToString(svg);
  var blob = new Blob([svgString], {type: 'image/svg+xml;charset=utf-8'});
  fileDownload(blob, title+'.svg');
  window.dispatchEvent(new Event('resize'));
}

function fileDownload(blob,name){
  if(window.navigator.msSaveBlob){
    window.navigator.msSaveBlob(blob, name);
  }else{
    var reader = new FileReader();
    reader.onload = function (event) {
      var save = document.createElement('a');
      save.href = event.target.result;
      save.target = '_blank';
      save.download = name;
      var clicEvent = new MouseEvent('click', {
        'view': window,
        'bubbles': true,
        'cancelable': true
      });
      save.dispatchEvent(clicEvent);
      (window.URL || window.webkitURL).revokeObjectURL(save.href);
    };
    reader.readAsDataURL(blob);
  }
}

function buttonLeftPos(sel){
  return (parseInt(sel.select("svg").style("width"))-20) + "px";
}

function iconButton(sel,alt,src,title,job,style){
    sel.append("img")
      .attr("class","icon")
      .attr("alt", alt)
      .style("width", "14px")
      .style("height", "14px")
      .attr("src", src)
      .attr("title", title)
      .style("cursor","pointer")
      .style(style)
      .on("click", job);
}

function formatter(d){
  if(typeof d == 'number'){
    var dabs = Math.abs(d);
    if((dabs>0 && dabs<1e-3) || dabs>1e+5)
      d = d.toExponential(3);
    else
      d = Math.round(d * 1000) / 1000;
  }
  return d;
}

(function() {

// Inspired by http://informationandvisualization.de/blog/box-plot
d3.box = function() {
  var width = 1,
      height = 1,
      duration = 0,
      domain = null,
      value = Number,
      whiskers = boxWhiskers,
      quartiles = boxQuartiles,
      tickFormat = null;

  // For each small multiple…
  function box(g) {
    g.each(function(d, i) {
      d = d.map(value).sort(d3.ascending);
      var g = d3.select(this),
          n = d.length,
          min = d[0],
          max = d[n - 1];

      // Compute quartiles. Must return exactly 3 elements.
      var quartileData = d.quartiles = quartiles(d);

      // Compute whiskers. Must return exactly 2 elements, or null.
      var whiskerIndices = whiskers && whiskers.call(this, d, i),
          whiskerData = whiskerIndices && whiskerIndices.map(function(i) { return d[i]; });

      // Compute outliers. If no whiskers are specified, all data are "outliers".
      // We compute the outliers as indices, so that we can join across transitions!
      var outlierIndices = whiskerIndices
          ? d3.range(0, whiskerIndices[0]).concat(d3.range(whiskerIndices[1] + 1, n))
          : d3.range(n);

      // Compute the new x-scale.
      var x1 = d3.scale.linear()
          .domain(domain && domain.call(this, d, i) || [min, max])
          .range([height, 0]);

      // Retrieve the old x-scale, if this is an update.
      var x0 = this.__chart__ || d3.scale.linear()
          .domain([0, Infinity])
          .range(x1.range());

      // Stash the new scale.
      this.__chart__ = x1;

      // Note: the box, median, and box tick elements are fixed in number,
      // so we only have to handle enter and update. In contrast, the outliers
      // and other elements are variable, so we need to exit them! Variable
      // elements also fade in and out.

      // Update center line: the vertical line spanning the whiskers.
      var center = g.selectAll("line.center")
          .data(whiskerData ? [whiskerData] : []);

      center.enter().insert("line", "rect")
          .attr("class", "center")
          .attr("x1", width / 2)
          .attr("y1", function(d) { return x0(d[0]); })
          .attr("x2", width / 2)
          .attr("y2", function(d) { return x0(d[1]); })
          .style("opacity", 1e-6)
        .transition()
          .duration(duration)
          .style("opacity", 1)
          .attr("y1", function(d) { return x1(d[0]); })
          .attr("y2", function(d) { return x1(d[1]); });

      center.transition()
          .duration(duration)
          .style("opacity", 1)
          .attr("y1", function(d) { return x1(d[0]); })
          .attr("y2", function(d) { return x1(d[1]); });

      center.exit().transition()
          .duration(duration)
          .style("opacity", 1e-6)
          .attr("y1", function(d) { return x1(d[0]); })
          .attr("y2", function(d) { return x1(d[1]); })
          .remove();

      // Update innerquartile box.
      var box = g.selectAll("rect.boxplot")
          .data([quartileData]);

      box.enter().append("rect")
          .attr("class", "boxplot")
          .attr("x", 0)
          .attr("y", function(d) { return x0(d[2]); })
          .attr("width", width)
          .attr("height", function(d) { return x0(d[0]) - x0(d[2]); })
        .transition()
          .duration(duration)
          .attr("y", function(d) { return x1(d[2]); })
          .attr("height", function(d) { return x1(d[0]) - x1(d[2]); });

      box.transition()
          .duration(duration)
          .attr("y", function(d) { return x1(d[2]); })
          .attr("height", function(d) { return x1(d[0]) - x1(d[2]); });

      // Update median line.
      var medianLine = g.selectAll("line.median")
          .data([quartileData[1]]);

      medianLine.enter().append("line")
          .attr("class", "median")
          .attr("x1", 0)
          .attr("y1", x0)
          .attr("x2", width)
          .attr("y2", x0)
        .transition()
          .duration(duration)
          .attr("y1", x1)
          .attr("y2", x1);

      medianLine.transition()
          .duration(duration)
          .attr("y1", x1)
          .attr("y2", x1);

      // Update whiskers.
      var whisker = g.selectAll("line.whisker")
          .data(whiskerData || []);

      whisker.enter().insert("line", "circle, text")
          .attr("class", "whisker")
          .attr("x1", 0)
          .attr("y1", x0)
          .attr("x2", width)
          .attr("y2", x0)
          .style("opacity", 1e-6)
        .transition()
          .duration(duration)
          .attr("y1", x1)
          .attr("y2", x1)
          .style("opacity", 1);

      whisker.transition()
          .duration(duration)
          .attr("y1", x1)
          .attr("y2", x1)
          .style("opacity", 1);

      whisker.exit().transition()
          .duration(duration)
          .attr("y1", x1)
          .attr("y2", x1)
          .style("opacity", 1e-6)
          .remove();

      // Update outliers.
      var outlier = g.selectAll("circle.outlier")
          .data(outlierIndices, Number);

      outlier.enter().insert("circle", "text")
          .attr("class", "outlier")
          .attr("r", 5)
          .attr("cx", width / 2)
          .attr("cy", function(i) { return x0(d[i]); })
          .style("opacity", 1e-6)
        .transition()
          .duration(duration)
          .attr("cy", function(i) { return x1(d[i]); })
          .style("opacity", 1);

      outlier.transition()
          .duration(duration)
          .attr("cy", function(i) { return x1(d[i]); })
          .style("opacity", 1);

      outlier.exit().transition()
          .duration(duration)
          .attr("cy", function(i) { return x1(d[i]); })
          .style("opacity", 1e-6)
          .remove();

      // Compute the tick format.
      var format = tickFormat || x1.tickFormat(8);

      // Update box ticks.
      var boxTick = g.selectAll("text.boxplot")
          .data(quartileData);

      boxTick.enter().append("text")
          .attr("class", "boxplot")
          .attr("dy", ".3em")
          .attr("dx", function(d, i) { return i & 1 ? 6 : -6 })
          .attr("x", function(d, i) { return i & 1 ? width : 0 })
          .attr("y", x0)
          .attr("text-anchor", function(d, i) { return i & 1 ? "start" : "end"; })
          .text(format)
        .transition()
          .duration(duration)
          .attr("y", x1);

      boxTick.transition()
          .duration(duration)
          .text(format)
          .attr("y", x1);

      // Update whisker ticks. These are handled separately from the box
      // ticks because they may or may not exist, and we want don't want
      // to join box ticks pre-transition with whisker ticks post-.
      var whiskerTick = g.selectAll("text.whisker")
          .data(whiskerData || []);

      whiskerTick.enter().append("text")
          .attr("class", "whisker")
          .attr("dy", ".3em")
          .attr("dx", 6)
          .attr("x", width)
          .attr("y", x0)
          .text(format)
          .style("opacity", 1e-6)
        .transition()
          .duration(duration)
          .attr("y", x1)
          .style("opacity", 1);

      whiskerTick.transition()
          .duration(duration)
          .text(format)
          .attr("y", x1)
          .style("opacity", 1);

      whiskerTick.exit().transition()
          .duration(duration)
          .attr("y", x1)
          .style("opacity", 1e-6)
          .remove();
    });
    d3.timer.flush();
  }

  box.width = function(x) {
    if (!arguments.length) return width;
    width = x;
    return box;
  };

  box.height = function(x) {
    if (!arguments.length) return height;
    height = x;
    return box;
  };

  box.tickFormat = function(x) {
    if (!arguments.length) return tickFormat;
    tickFormat = x;
    return box;
  };

  box.duration = function(x) {
    if (!arguments.length) return duration;
    duration = x;
    return box;
  };

  box.domain = function(x) {
    if (!arguments.length) return domain;
    domain = x == null ? x : d3.functor(x);
    return box;
  };

  box.value = function(x) {
    if (!arguments.length) return value;
    value = x;
    return box;
  };

  box.whiskers = function(x) {
    if (!arguments.length) return whiskers;
    whiskers = x;
    return box;
  };

  box.quartiles = function(x) {
    if (!arguments.length) return quartiles;
    quartiles = x;
    return box;
  };

  return box;
};

function boxWhiskers(d) {
  return [0, d.length - 1];
}

function boxQuartiles(d) {
  return [
    d3.quantile(d, .25),
    d3.quantile(d, .5),
    d3.quantile(d, .75)
  ];
}

})();
