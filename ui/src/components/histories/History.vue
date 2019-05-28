<template>
    <div class="galaxy-history">
        <span class="galaxy-history-label" v-bind:contenteditable="editing_name">{{ model.name }}</span>
        <span class="galaxy-history-state">{{ model.state }}</span>
        <time class="galaxy-history-updated" v-bind:datetime="model.update_time">{{ (new Date(model.update_time)).toLocaleDateString(undefined, { weekday: 'long', year: 'numeric', month: 'long', day: 'numeric', hour: 'numeric', minute: 'numeric', second: 'numeric' }) }}</time>
        <slot v-bind:model="model" class="galaxy_history_slot"></slot>
        <div class="galaxy-history-functions">
            <slot name="functions" v-bind="self" />
            <!--TODOa @click.stop.prevent="download" href="">Prepare Download</a-->
            <a @click.stop.prevent="remove" href="">Remove</a>
        </div>
        <slot name="contents" v-bind="self" class="galaxy-history-contents"/>
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

        },
        data() {return {
            self: this,
            editing_name: false,
        }},
        methods: {
            remove() {
                this.model.delete({query: {purge: true}}); //TODO Remove purge=True?
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