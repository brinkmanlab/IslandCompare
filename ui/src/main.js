import Vue from 'vue'
import AsyncComputed from 'vue-async-computed'
import App from './App.vue'
import store from './store'

//import i18n from './i18n'

Vue.config.productionTip = false;

Vue.use(AsyncComputed);

new Vue({
  store,
//  i18n,
  render: h => h(App),
}).$mount('#app');
