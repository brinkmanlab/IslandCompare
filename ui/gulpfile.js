const { src, dest } = require('gulp');
const md = require('gulp-markdownit');
const rename = require('gulp-rename');
const header = require('gulp-header');

const MD_PLUGINS = ['markdown-it-deflist'];

function defaultTask() {
    // Render static markdown to html fragments
    return src('static/*.md')
        .pipe(md({/*options: {html: true, linkify: true}*/ plugins: MD_PLUGINS}))
        .pipe(rename({
            extname: ".htm",
        }))
        .pipe(header('<!-- DO NOT MODIFY! THIS CONTENT WAS DYNAMICALLY GENERATED AND ANY CHANGES WILL BE OVERWRITTEN. See ../static/*.md -->\n'))
        .pipe(dest('public/'));
}

exports.default = defaultTask;
