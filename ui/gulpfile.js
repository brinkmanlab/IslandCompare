const { src, dest, parallel } = require('gulp');
const md = require('gulp-markdownit');
const rename = require('gulp-rename');
const header = require('gulp-header');
const markdownToJSON = require('gulp-markdown-to-json');
//const ListStream = require('list-stream');
const concatJson = require('gulp-concat-json-to-array');

const MD_PLUGINS = ['markdown-it-deflist'];
const MD_CONFIG = {options: {html: true/*, linkify: true*/}, plugins: MD_PLUGINS};

const md_inst = require('markdown-it')(MD_CONFIG.preset || 'default', MD_CONFIG.options || {});

MD_PLUGINS.forEach(plugin => {
    md_inst.use(require(plugin), md_inst.options);
});

function pages() {
    // Render static markdown to html fragments
    return src('static/*.md')
        .pipe(md(MD_CONFIG))
        .pipe(rename({
            extname: ".htm",
        }))
        .pipe(header('<!-- DO NOT MODIFY! THIS CONTENT WAS DYNAMICALLY GENERATED AND ANY CHANGES WILL BE OVERWRITTEN. See ../static/*.md -->\n'))
        .pipe(dest('src/assets/'));
}

function news() {
    // Render all news markdown to json, sorted by date
    return src('static/news/*.md')
        //.pipe(ListStream.obj())
        .pipe(markdownToJSON(md_inst.render.bind(md_inst), 'news.json'))
        .pipe(concatJson('news.json', (data)=>new Buffer(JSON.stringify(data.map(item=>{item.date = new Date(item.date); return item}).sort((a,b)=>b.date-a.date))))) // TODO list-streams broken for some reason, subbing this in
        .pipe(dest('src/assets/'));
}

function slides() {
    // Render all slides to json, sorted by file name
    return src('static/slides/*.md')
    //.pipe(ListStream.obj())
        .pipe(markdownToJSON(md_inst.render.bind(md_inst), 'slides.json'))
        .pipe(concatJson('slides.json', (data)=>new Buffer(JSON.stringify(data)))) // TODO list-streams broken for some reason, subbing this in
        .pipe(dest('src/assets/'));
}

exports.default = parallel(pages, news, slides);
