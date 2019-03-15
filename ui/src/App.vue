<template>
  <div id="app">
    <FileUploader/>
  </div>
</template>

<script>
import FileUploader from './components/JobRunner.vue'

export default {
  name: 'app',
  components: {
    FileUploader
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
