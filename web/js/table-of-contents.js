var ToC = "<nav role='navigation' class='table-of-contents'>" + 
			"<span class=\"font-21\">Table of Contents:</span>" +
   	 		"<ul>";

var newLine, el, title, link;

$("span.font-21").each(function() {

  el = $(this);
  title = el.text();
  link = "#" + el.attr("id");

  newLine =
    "<li><span class=\"font\">" +
      "<a href='" + link + "'>" +
        title +
      "</a>" +
    "</span></li>";
  ToC += newLine;

});

ToC +=
   "</ul>" +
  "</nav>";
	
$(".cont-52").prepend(ToC);