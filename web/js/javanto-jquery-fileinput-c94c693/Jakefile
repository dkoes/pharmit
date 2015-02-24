desc("This is the default task.");
task("default", function (params) {
  var fs = require("node-fs");
  var dot = require("dot");
  var uglify = require("uglify-js");

  package(function(pkg) {
    console.log("Building " + pkg.title + " " + pkg.version);
    fs.readFile("./jquery.fileinput.js", "utf8", function(err, data) {
      if (err) throw err;
      var cmp = dot.template(data, {evaluate: /\{\{([\s\S]+?)\}\}/g,
        interpolate: /\{\{=([\s\S]+?)\}\}/g,
        encode: /\{\{!([\s\S]+?)\}\}/g,
        use: /\{\{#([\s\S]+?)\}\}/g, //compile time evaluation
        define: /\{\{##\s*([\w\.$]+)\s*(\:|=)([\s\S]+?)#\}\}/g, //compile time defs
        strip: false,
        varname: "it"})(pkg);
      var ast = uglify.parse(cmp);
      ast.figure_out_scope();
      ast = ast.transform(uglify.Compressor());
      // need to figure out scope again so mangler works optimally
      ast.figure_out_scope();
      ast.compute_char_frequency();
      ast.mangle_names();

      var min = ast.print_to_string({ comments: true });
      var buildDirectory = pkg.build.directory? dot.template(pkg.build.directory)(pkg) : "dist";
      try {
        fs.rmdirSync("./" + buildDirectory);
      } catch (e) {
      }
      fs.mkdirSync("./" + buildDirectory, 0777, true);
      console.log("Writing to " + buildDirectory);
      var buildFinalName = pkg.build.final_name ? dot.template(pkg.build.final_name)(pkg) : pkg.name + "-" + pkg.version;
      fs.writeFileSync("./" + buildDirectory + "/" + buildFinalName + ".min.js", min);
      fs.writeFileSync("./" + buildDirectory + "/" + buildFinalName + ".js", cmp);
      complete();
    });
  });
}, true);

task("publish", ["default"], function (params) {
  var gits = require("gits");

  package(function(pkg) {
    var tagName = pkg.version;
    gits.git(".", ["tag", "-a", "-m", "Publishing new version", tagName], function() {
      gits.git(".", ["push", "upstream", "--tags"], function() {
        console.log("Created remote tag " + tagName);
      });
    });
  });
});

function package(callback) {
  var fs = require("node-fs");

  fs.readFile("./fileinput.jquery.json", "utf8", function(err, data) {
    if (err) throw err;
    callback(JSON.parse(data));
  });
}
