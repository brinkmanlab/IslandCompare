# Contributing

Contributions are welcome.

Submit contributions as a pull request on GitHub.

## Editing Static Content
All static content is placed in `ui/static` as markdown formatted files.

### News
News items can be added to the home page by creating new markdown files in `ui/static/news/`.
The file names should be of the form `{YEAR}-{MONTH}-{DAY}-{SLUG}.md`. 
 - YEAR - 4 digit year (2019)
 - MONTH - 2 digit month (06)
 - DAY - 2 digit day (02)
 - SLUG - a short description of the content, words separated by dashes (`-`). The slug must be unique for all news items.

Every file must contain a header referred to as the 'frontmatter'.
The frontmatter is of the form:

```markdown

---
title: Title to be displayed above the body
date: News date, formatted the same as the file name.
slug: The same slug used in the file name
---

News body content goes here.

``` 

News items are sorted by the date specified in the frontmatter.

See [this markdown cheatsheet](https://github.com/adam-p/markdown-here/wiki/Markdown-Cheatsheet) for more information on formatting the body text.


### Slides
Slides can be added to the home page carousel by creating new markdown files in `ui/static/slides/`.
Each slide has a title, body of text, and accompanying image.
Slide images are stored in `ui/static/slides/images/`.

Slide filenames are not important other than they must end in `.md`, not contain spaces, and are ordered alphanumerically. 
It is recommended that all file names are prefixed with a number so that order is explicitly clear.

Slide markdown files must contain frontmatter of the form:

```markdown

---
title: Title to be displayed above the body
slug: a short description of the content, words separated by dashes (`-`). The slug must be unique for all slides.
img: filename of the image in `ui/static/slides/images/` 
alt: a short description of the image for the blind
---

Slide body text goes here

```

See [this markdown cheatsheet](https://github.com/adam-p/markdown-here/wiki/Markdown-Cheatsheet) for more information on formatting the body text.