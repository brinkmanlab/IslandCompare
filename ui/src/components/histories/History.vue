<template>
    <div class="History">
        <span v-bind:class="formatter('galaxy-history-label', model.name, self)">{{ model.name }}</span>
        <span v-bind:class="formatter('galaxy-history-state', model.state, self)">{{ model.state }}</span>
        <time v-bind:class="formatter('galaxy-history-updated', model.update_time, self)" v-bind:datetime="model.update_time">{{ (new Date(model.update_time)).toLocaleDateString(undefined, { weekday: 'long', year: 'numeric', month: 'long', day: 'numeric', hour: 'numeric', minute: 'numeric', second: 'numeric' }) }}</time>
        <slot v-bind:model="model" v-bind:class="formatter('galaxy_history_slot', model, self)"></slot>
        <div v-bind:class="formatter('galaxy-history-functions', self, self)">
            <slot name="functions" v-bind="self" />
            <!--TODOa @click.stop.prevent="download" href="">Prepare Download</a-->
            <a @click.stop.prevent="remove" href="">Remove</a>
        </div>
        <slot name="contents" v-bind="self" v-bind:class="formatter('galaxy-history-contents')"/>
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
            },
            formatter: {
                type: Function,
                default: c=>c,
            }
        },
        data() {return {
            self: this,
        }},
        methods: {
            remove() {
                galaxy.histories.History.$delete({params: {id: this.model.id}}); //TODO Add purge=True?
            }
        },
        mounted() {
            if (!this.model.constructor.end_states.includes(this.model.state)) {
                this.model.start_polling(()=>{
                    if (this.model.constructor.end_states.includes(this.model.state)) {
                        this.$emit('history-completed', this);
                        return true;
                    }
                    return false;
                });
            }
        },
        beforeDestroy() {
            if (this.model) this.model.stop_polling();
        }
    }
</script>

<style scoped>

</style>