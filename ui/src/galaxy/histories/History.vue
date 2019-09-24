<template xmlns:v-slot="http://www.w3.org/1999/XSL/Transform">
    <div class="galaxy-history">
        <EditableLabel class="galaxy-history-label" @update="update_label" v-bind:value="model.name" placeholder="Enter a label" ref="label"></EditableLabel>
        <span class="galaxy-history-state">{{ model.state }}</span>
        <time class="galaxy-history-updated" v-bind:datetime="model.update_time">{{ (new Date(model.update_time)).toLocaleDateString(undefined, { weekday: 'long', year: 'numeric', month: 'long', day: 'numeric', hour: 'numeric', minute: 'numeric', second: 'numeric' }) }}</time>
        <slot v-bind:model="model" class="galaxy_history_slot"></slot>
        <HistoryFunctions v-bind:item="self">
            <template v-slot:default="slot">
                <slot name="functions" v-bind="self" />
                <!--TODOa @click.stop.prevent="download" href="">Prepare Download</a-->
                <RenameHistory v-bind:item="slot.item" v-on:galaxy-history-rename="$refs.label.start_edit()"/>
                <RemoveHistory v-bind:item="slot.item" v-on="$listeners"/>
            </template>
        </HistoryFunctions>
        <slot name="contents" v-bind="self" class="galaxy-history-contents"/>
    </div>
</template>

<script>
    import { History } from "@/galaxy/api/histories";
    import EditableLabel from "@/galaxy/misc/EditableLabel";
    import HistoryFunctions from "@/galaxy/histories/HistoryFunctions"
    import RemoveHistory from "@/galaxy/histories/HistoryFunctions/Remove";
    import RenameHistory from "@/galaxy/histories/HistoryFunctions/Rename";
    export default {
        name: "History",
        components: {HistoryFunctions, RemoveHistory, RenameHistory, EditableLabel},
        props: {
            model: {
                type: History,
                required: true,
            },
        },
        data() {return {
            self: this,
        }},
        methods: {
            update_label(value) {
                this.model.name = value;
                this.model.post(['name']);
            },
        },
        mounted() {
            const self = this;
            this.model.poll_state(()=>{self.$emit('history-completed', self); return true;});
        },
        beforeDestroy() {
            if (this.model) this.model.stop_polling();
        }
    }
</script>

<style scoped>

</style>