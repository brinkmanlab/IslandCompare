<template xmlns:v-slot="http://www.w3.org/1999/XSL/Transform">
    <b-container class="galaxy-history-item">
        <slot name="before"/>
        <span class="galaxy-history-item-hid">{{ model.hid }}</span>
        <EditableLabel class="galaxy-history-item-name" @update="update_label" v-bind:value="model.name" placeholder="Enter a name to identify this dataset" ref="label"></EditableLabel>
        <!-- TODO add more progress states depending on hda state -->
        <b-progress v-bind:class="'galaxy-history-item-progress w-100 '+this.model.state"
                    v-if="model.upload_progress<100"
                    v-bind:max="100"
                    variant="info" striped animated>
            <b-progress-bar class="galaxy-history-item-progressbar" v-bind:value="model.upload_progress"></b-progress-bar>
        </b-progress>
        <slot></slot>
        <HistoryItemFunctions v-bind:item="this" v-on:galaxy-history-item-rename="$refs.label.start_edit()" v-on="$listeners">
            <template v-slot:default="slot">
                <slot name="history_item_functions" v-bind:item="slot.item" v-bind:listeners="$listeners"></slot>
            </template>
        </HistoryItemFunctions>
        <slot name="after"/>
    </b-container>
</template>

<script>
    import {HistoryDatasetAssociation, HistoryDatasetCollectionAssociation} from "@/galaxy/api/history_contents";
    import HistoryItemFunctions from "./HistoryItemFunctions";
    import EditableLabel from "@/galaxy/misc/EditableLabel";

    export default {
        name: "HistoryItem",
        components: {
            EditableLabel,
            HistoryItemFunctions,
        },
        data: ()=>{return{
        }},
        props: {
            model: {
                type: [HistoryDatasetAssociation, HistoryDatasetCollectionAssociation],
                required: true,
            },
        },
        computed: {
        },
        methods: {
            update_label(evt, value) {
                this.model.name = value;
                this.model.upload(['name']);
            }
        },
        mounted() {
        },
    }
</script>

<style scoped>
    .galaxy-history-item {
        display: flex;
        flex-direction: row;
        align-items: center;
        width: 100%;
    }

    .galaxy-history-item >>> .galaxy-history-item-functions {
        height: 1.1em;
    }

</style>