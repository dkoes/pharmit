# jQuery Fileinput Plugin

Normally, you can't use CSS to style HTML file inputs. They'll always look [fugly](http://www.urbandictionary.com/define.php?term=fugly). With jQuery Fileinput you can!

By wrapping the file input on a div, adding a replacement element to provide the visuals and then making the original input invisible. Clicks to the replacement element are propagated to the original input on the browsers that support it, otherwise a transparent file input is placed on top of the replacement.

## Downloads
* [Minified](https://raw.github.com/javanto/jquery-fileinput/3.1.0/dist/jquery.fileinput.min.js)
* [Full version](https://raw.github.com/javanto/jquery-fileinput/3.1.0/dist/jquery.fileinput.min.js)

## Examples

### Trivial
    ...
    <script type="text/javascript" src="https://ajax.googleapis.com/ajax/libs/jquery/1.9.0/jquery.min.js"></script>
    <script type="text/javascript" src="https://raw.github.com/javanto/jquery-fileinput/3.1.0/dist/jquery.fileinput.min.js"></script>
    <script type="text/javascript">
    //<![CDATA[
    $().ready(function() {
      $("input[type=file]").fileinput();
    });
    //]]>
    </script>
    ...
    <input type="file" name="foo" />

### With replacement HTML
    ...
    <script type="text/javascript" src="https://ajax.googleapis.com/ajax/libs/jquery/1.9.0/jquery.min.js"></script>
    <script type="text/javascript" src="https://raw.github.com/javanto/jquery-fileinput/3.1.0/dist/jquery.fileinput.min.js"></script>
    <script type="text/javascript">
    //<![CDATA[
    $().ready(function() {
      $("input[type=file]").fileinput("<button>Browse...</button>");
    });
    //]]>
    </script>
    ...
    <input type="file" name="foo" />

### With replacement jQuery selector
    ...
    <script type="text/javascript" src="https://ajax.googleapis.com/ajax/libs/jquery/1.9.0/jquery.min.js"></script>
    <script type="text/javascript" src="https://raw.github.com/javanto/jquery-fileinput/3.1.0/dist/jquery.fileinput.min.js"></script>
    <script type="text/javascript">
    //<![CDATA[
    $().ready(function() {
      $("input[type=file]").fileinput($("button").remove());
    });
    //]]>
    </script>
    ...
    <input type="file" name="foo" />
    <button>Browse...</button>

## API

### .fileinput([replacement])

* **replacement** Plain HTML, or a jQuery selector for the element intended to replace the file input.

### Notes

With **Opera** (and **Firefox** versions before 4) we can't support the [native CSS pseudo classes](http://www.w3schools.com/css/css_pseudo_classes.asp). Instead we'll provide them equivalent "pseudo pseudo classes" which are actually just normal CSS classes.

* **.hover** for :hover
* **.focus** for :focus
* **.active** for :active

## Demo

At [jsFiddle](http://jsfiddle.net/hleinone/UF4nr/).

## Requirements

* [jQuery](http://jquery.com/)
 * 1.4+

## Brower support

### Tested on

* Chrome
 * 24.0.1312.57 m
* Chromium
 * 24.0.1312.56
* Firefox
 * 18.0.2
* Internet Explorer
 * 10.0.9200.16466
* Opera
 * 12.14
* Safari
 * 5.1.7 (7534.57.2)

## Developing

1. Fork
1. Make your modifications
1. Make sure you have [node](http://nodejs.org/) and [npm](http://npmjs.org/) installed
1. Install the dependencies required for the build: `npm install node-fs dot uglify-js gits`
1. Install [jake](https://github.com/mde/jake): `npm install jake -g`
1. Compile and minify: `jake` 
1. Test
1. Commit and create a pull request
1. ???
1. Profit!
