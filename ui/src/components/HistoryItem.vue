<template>
    <li>
        <slot name="before"/>
        <span class="name">{{hid}}: {{name}}</span>
        <progress max="100" v-bind:value="progress" v-if="progress>0"></progress>
        <slot></slot>
        <HistoryItemFunctions v-bind:item="this">
            <template v-slot:default="slot">
                <slot name="history_item_functions" v-bind:item="slot.item"></slot>
                <RemoveHistoryItem v-bind:item="slot.item" />
            </template>
        </HistoryItemFunctions>
        <slot name="after"/>
    </li>
</template>

<script>
    import axios from 'axios';
    import HistoryItemFunctions from "./HistoryItemFunctions";
    import RemoveHistoryItem from "./HistoryItemFunctions/Remove";

    export default {
        name: "HistoryItem",
        components: {
            HistoryItemFunctions,
            RemoveHistoryItem,
        },
        data: ()=>{return{
            hid: 0,
            name: '',
            progress: 0,
        }},
        props: {
            id: {
                type: String,
                required: true,
            },
            initial_data: {
                type: Object,
            },

        },
        computed: {
        },
        methods: {
            refresh() {
                axios.get(`/api/histories/${this.history_id}/contents/`, {
                    params: {
                        key: this.$store.state.galaxy.api_key,
                    }

                }).then(response => {
                    Object.assign(this.$data, response.data);
                }).catch(error => {
                    //TODO
                    console.log(error.status); // eslint-disable-line no-console
                });
            }
        },
        mounted() {
            if (this.$props.initial_data)
                Object.assign(this.$data, this.$props.initial_data);
            else {
                this.refresh();
            }
        },
    }
</script>

<style scoped>

</style>