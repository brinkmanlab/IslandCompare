<template>
  <div id="app">
    <FileUploader @toast="toast"/>
    <div ref="toast" class="toast"></div>
  </div>
</template>

<script>
import FileUploader from './components/JobRunner.vue'

export default {
  name: 'app',
  components: {
    FileUploader
  },
  methods: {
    toast(e) {
      //TODO
      console.log(e); // eslint-disable-line no-console
    },
  }
}
</script>

<style>
#app {
  font-family: 'Avenir', Helvetica, Arial, sans-serif;
  -webkit-font-smoothing: antialiased;
  -moz-osx-font-smoothing: grayscale;
}
</style>
