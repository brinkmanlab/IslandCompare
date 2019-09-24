<template>
    <b-container class="galaxy-reference-genomes">
        <b-row>
            <b-input-group size="sm" prepend="Filter">
                <b-form-input
                        v-bind:value="filter"
                        @input="update_filter"
                        type="search"
                        placeholder="Type to Search"
                />
                <b-input-group-append>
                    <b-button :disabled="!filter" @click="filter = ''">Clear</b-button>
                </b-input-group-append>
            </b-input-group>
        </b-row>
        <b-table hover small borderless selectable show-empty
                 select-mode="single"
                 v-bind:items="genomes"
                 v-bind:busy="isLoading"
                 v-bind:fields="[{key: 'name', sortable: true}]"
                 v-bind:current-page="currentPage"
                 v-bind:per-page="perPage"
                 v-bind:filter="filter"
                 v-bind:filter-function="filterFunc"
                 sort-by="name"
                 primary-key="id"
                 @row-selected="row_selected"
                 thead-class="hidden_header"
                 ref="table"
        >
            <template slot="table-busy">
                <div class="text-center">
                    <b-spinner class="align-middle"></b-spinner>
                    <strong>Loading...</strong>
                </div>
            </template>
        </b-table>
        <b-row align-h="center">
            <b-pagination
                    v-model="currentPage"
                    :per-page="perPage"
                    :total-rows="totalRows"
                    limit="8"
            ></b-pagination>
        </b-row>
    </b-container>
</template>

<script>
    import { Genome } from "@/galaxy/api/genomes";

    export default {
        name: "ReferenceGenomes",
        props: {
            perPage: {
                type: Number,
                default: 10,
            },
            value: {
                type: String,
                default: '',
            },
            filter_delay: {
                type: Number,
                default: 500,
            }
        },
        data(){return{
            currentPage: 1,
            filter: "",
            filter_delay_handle: null,
        }},
        methods: {
            row_selected(list) {
                if (list.length)
                    this.$emit('input', list[0]);
                else
                    this.$emit('input', '');
            },
            clearSelected() { this.$refs.table.clearSelected(); this.$emit('input', ''); },
            filterFunc(row, filter) {
                return row.name.includes(filter);
            },
            update_filter(filter) {
                // Delay the filter slightly to reduce lag and allow the user to continue typing
                if (this.filter_delay_handle !== null) window.clearTimeout(this.filter_delay_handle);
                this.filter_delay_handle = window.setTimeout(t=>t.filter = filter, this.filter_delay, this);
            }
        },
        computed: {
            totalRows() {
                return this.genomes.length;
            },
            isLoading() {
                return this.genomes.length === 0;
            }
        },
        asyncComputed: {
            genomes: {
                async get() {
                    let genomes = Genome.all();
                    if (genomes.length === 0) {
                        // Fetch once
                        await Genome.$fetch();
                        genomes = Genome.all();
                    }
                    return genomes;
                },
                default: [],
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
    }
</script>

<style>
    .galaxy-reference-genomes .hidden_header {
        display: none;
    }

    .galaxy-history-contents .row:first-child {
        flex-wrap: nowrap;
        align-items: center;
    }
</style>