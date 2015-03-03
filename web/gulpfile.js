var gulp = require('gulp'),
    jshint = require('gulp-jshint'),
    concat = require('gulp-concat'),
    uglify = require('gulp-uglify');
    rename = require('gulp-rename');
 
// Lint JS
gulp.task('lint', function() {
  gulp.src('./src/*.js')
    .pipe(jshint({'latedef':'nofunc'}))
    .pipe(jshint.reporter('default'));
});
 
// Concat & Minify JS
gulp.task('minify', function(){
    gulp.src('./src/*.js')
        .pipe(concat('pharmit.js'))
        .pipe(gulp.dest('js'))
        .pipe(rename('pharmit.min.js'))
        .pipe(uglify())
        .pipe(gulp.dest('js'));
});
 
// Default
gulp.task('default', ['lint','minify'] );

// Watch Our Files
gulp.task('watch', function() {
  gulp.watch('src/*.js', ['lint', 'minify']);
});
