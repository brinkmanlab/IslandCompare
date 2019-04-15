import Vue from 'vue'
import AsyncComputed from 'vue-async-computed'
import App from './App.vue'
import { getStore } from './store'

//TODO import i18n from './i18n'

Vue.config.productionTip = false;

Vue.use(AsyncComputed);

getStore().then(store=>{
    new Vue({
        store,
        //  i18n,
        render: h => h(App),
    }).$mount('#app');
});