var gulp = require('gulp'),
    jshint = require('gulp-jshint'),
    concat = require('gulp-concat'),
    uglify = require('gulp-uglify');
    rename = require('gulp-rename');
const { series } = require('gulp'); 

psrc = ['src/feature.js',
        'src/hex_md5.js',
        'src/jquery.fileinput.js',
        'src/minresults.js',
        'src/msg.js',
        'src/query.js',
        'src/results.js',
        'src/searchresults.js',
        'src/spectrum.js',
        'src/viewer.js',
        'src/main.js'];
        
// Lint JS
gulp.task('lint', function() {
  return gulp.src(psrc)
    .pipe(jshint({'latedef':'nofunc'}))
    .pipe(jshint.reporter('default'));
});
 
// Concat & Minify JS
gulp.task('minify', function(){
    return gulp.src(psrc)
        .pipe(concat('pharmit.js'))
        .pipe(gulp.dest('js'))
        .pipe(rename('pharmit.min.js'))
        .pipe(uglify())
        .pipe(gulp.dest('js'));
});
 
// Default
gulp.task('default', series('lint','minify') );

// Watch Our Files
gulp.task('watch', function() {
  gulp.watch('src/*.js', ['lint', 'minify']);
});
