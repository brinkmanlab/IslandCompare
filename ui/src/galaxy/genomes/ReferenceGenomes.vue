<template>
    <b-container fluid>
        <b-row>
            <b-col md="6">
                <b-form-group label-cols-sm="3" label="Filter" class="mb-0">
                    <b-input-group>
                        <b-form-input v-model="filter" placeholder="Type to Search"></b-form-input>
                        <b-input-group-append>
                            <b-button :disabled="!filter" @click="filter = ''">Clear</b-button>
                        </b-input-group-append>
                    </b-input-group>
                </b-form-group>
            </b-col>
        </b-row>
        <b-table hover small borderless selectable
                 select-mode="single"
                 v-bind:items="genomes"
                 v-bind:fields="[{key: 'name', sortable: true}]"
                 v-bind:current-page="currentPage"
                 v-bind:per-page="perPage"
                 v-bind:filter="filter"
                 primary-key="id"
                 @row-selected="row_selected"
        />
        <b-row>
            <b-col md="6" class="my-1">
                <b-pagination
                        v-model="currentPage"
                        :total-rows="totalRows"
                        :per-page="perPage"
                        class="my-0"
                ></b-pagination>
            </b-col>
        </b-row>
    </b-container>
</template>

<script>
    import { Genome } from "@/galaxy/api/genomes";
    export default {
        name: "ReferenceGenomes",
        props: {
            totalRows: {
                type: Number,
                default: 10,
            },
            perPage: {
                type: Number,
                default: 10,
            },
            value: {
                type: String,
                default: '',
            }
        },
        data(){return{
            currentPage: 0,
            filter: "",
        }},
        methods: {
            row_selected(list) {
                if (list.length)
                    this.$emit('input', list[0]);
                else
                    this.$emit('input', '');
            },
            clearSelected() { this.$refs.table.clearSelected(); this.$emit('input', ''); },
        },
        computed: {
            genomes() {
                return Genome.getAll();
            }
        },
        mounted() {
            if (this.value) {
                for (const [index, genome] of this.genomes().elements()) {
                    if (this.value === genome.id) {
                        this.$refs.table.selectRow(index);
                        break;
                    }
                }
            }
        },
        beforeMount() {
            Genome.$fetch()
        },
    }
</script>

<style scoped>

</style>