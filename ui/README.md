# IslandCompare

## Notes
- Datasets are stored in a history tagged 'user_data'
- Each job gets its own history tagged with the workflow id
- Be sure to update vue.config.js/publicPath if the web path changes
- Static content is stored in markdown files in ./static. See https://github.com/markdown-it/markdown-it for syntax extensions.
- See src/app.config.js for configuration parameters

## Directory
- src/components - items reusable for other galaxy projects
- src/IslandCompare - components specific to brinkman lab
- src/api - all models relating to the galaxy api. Many are stubbed and need to be filled out.
- src/galaxy.js - main import for galaxy api models and initialization
- src/auth.js - all code related to authenticating with the backend
- src/store.js - Vuex init code 
- src/main.js - main entry point of app
- src/app.js - all code specifically for loading and invoking IslandCompare state 
- gulpfile.js - responsible for all preprocessing before webpack is invoked


## Project setup
```
npm install
```

### Compiles and hot-reloads for development
```
npm run serve
```

### Compiles and minifies for production
```
npm run build
```

### Run your tests
```
npm run test
```

### Lints and fixes files
```
npm run lint
```

### Customize configuration
See [Configuration Reference](https://cli.vuejs.org/config/).
