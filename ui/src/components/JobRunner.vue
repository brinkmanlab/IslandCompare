<template>
    <div class="JobRunner">
        <HistoryContents v-if="model" v-bind:model="model" filter="dataset" @operation-fail="e=>this.$emit('toast', e)"/>
        <input type="submit" />
    </div>
</template>

<script>
    //import { Model as UserModel } from '@/api/users';
    import { Model as HistoryModel } from '@/api/histories'
    import HistoryContents from './HistoryContents';
    export default {
        name: "JobRunner",
        components: {
            HistoryContents
        },
        data: ()=>{return {
            model: null,
        }},
        computed: {
        },
        methods: {
        },
        mounted() {
            const self = this;
            HistoryModel.$fetch();
            this.model = HistoryModel.query().where('tags', tags=>tags.includes('user_data')).first();
            if (!this.model) {
                HistoryModel.$create({
                    data: {
                        name: "Uploaded data",
                    }
                }).then(response => {
                    self.model = HistoryModel.find(response.id);
                    self.model.tags.push('user_data');
                    self.model.$update({
                        params: {
                            id: self.model.id,
                        }
                    }).then((response)=>{
                        console.log("Updated"); //eslint-disable-line
                        console.log(response);//eslint-disable-line
                    }).catch((response)=>{
                        console.log("Failed");//eslint-disable-line
                        console.log(response);//eslint-disable-line
                    });
                });
            }
        },
    }
</script>

<style scoped>
    .JobRunner input[type=submit] {

    }
</style>