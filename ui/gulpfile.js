const { src, dest } = require('gulp');
const md = require('gulp-markdown-it');
const rename = require("gulp-rename");

function defaultTask() {
    // Render static markdown to html fragments
    return src('*.md')
        .pipe(md(/*{options: {html: true, linkify: true}}*/))
        .pipe(rename({
            extname: ".htm",
        }))
        .pipe(dest('../public/'));
}

exports.default = defaultTask
