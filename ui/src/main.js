import Vue from 'vue'
import App from './App.vue'
import store from './store'
import VueCookies from 'vue-cookies'

Vue.use(VueCookies);

// set default config
VueCookies.config('7d');

//import i18n from './i18n'

Vue.config.productionTip = false;

new Vue({
  store,
//  i18n,
  render: h => h(App),
}).$mount('#app');
