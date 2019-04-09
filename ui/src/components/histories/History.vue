<template>
    <div class="History">
        <span>{{ model.name }}</span>
        <span>{{ model.state }}</span>
        <time v-bind:datetime="model.update_time">{{ (new Date(model.update_time)).toLocaleDateString(undefined, { weekday: 'long', year: 'numeric', month: 'long', day: 'numeric', hour: 'numeric', minute: 'numeric', second: 'numeric' }) }}</time>
        <slot v-bind:model="model"></slot>
        <div class="functions">
            <slot name="functions" v-bind="this" />
            <!--TODOa @click.stop.prevent="download" href="">Prepare Download</a-->
            <a @click.stop.prevent="remove" href="">Remove</a>
        </div>
    </div>
</template>

<script>
    import * as galaxy from "@/galaxy";
    export default {
        name: "History",
        props: {
            model: {
                type: galaxy.histories.History,
                required: true,
            }
        },
        data() {return {
            pollHandle: null,
            stopPoll: ['ok','error'],
        }},
        methods: {
            remove() {
                galaxy.histories.History.$delete({params: {id: this.model.id}}); //TODO Add purge=True?
            }
        },
        mounted() {
            if (!this.stopPoll.includes(this.model.state)) {
                this.pollHandle = setInterval(()=>{
                    galaxy.histories.History.$get({params:{id: this.model.id}}).then(()=>{
                        if (this.stopPoll.includes(this.model.state)) {
                            clearInterval(this.pollHandle);
                            this.pollHandle = null;
                            this.$emit('history-completed', this);
                        }
                    });
                }, 10000);
            }
        },
        beforeDestroy() {
            if (this.pollHandle !== null) clearInterval(this.pollHandle);
        }
    }
</script>

<style scoped>

</style>