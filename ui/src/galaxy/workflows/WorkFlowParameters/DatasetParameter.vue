<template>
    <b-card class="galaxy-workflow-parameter-dataset" v-bind:header="label+(optional?' (Optional)': '')" v-bind:border-variant="validation_message ? 'danger' : 'default'">
        <HistoryContents ref="history_contents"
                         v-bind:historyPromise="historyPromise"
                         v-bind:value="value"
                         v-bind:filter="datasetFilter"
                         v-bind:accepted_upload_types="format"
                         @input="onInput"
                         @upload="upload=>this.$emit('upload', upload)"
        />
        <span v-if="mapped">
            <i class="icon-batch-mode"></i>
            This input expects a single selection. Separate jobs will be triggered for each dataset selection.
        </span>
        <b-card-footer v-if="validation_message" footer-text-variant="danger">
            <em>{{validation_message}}</em>
        </b-card-footer>
        <b-card-footer v-else footer-text-variant="info">
            {{ annotation }}
        </b-card-footer>
    </b-card>
</template>

<script>
    import HistoryContents from "../../histories/HistoryContents";
    export default {
        name: "DatasetParameter",
        components: { HistoryContents },
        props: {
            type: {
                type: String,
                required: true,
            },
            label: {
                type: String,
                required: true,
            },
            annotation: {
                type: String,
                default: '',
            },
            historyPromise: {
                type: Promise,
                required: true,
            },
            value: {
                type: Array,
                default(){return []},
            },
            format: {
                type: Array,
                default() { return [] },
            },
            optional: {
                type: Boolean,
                default: false,
            },
            validator: {
                type: Function,
                default(selection) {
                    return selection ? '' : 'Please select at least one dataset';
                }
            }
        },
        data() {return{
            selection: null,
            selected_models: [],
            validation_message: '',
        }},
        computed: {
            mapped() {
                return this.type === 'data_input' && this.selection && this.selection.src === 'new_collection';
            },
        },
        methods: {
            datasetFilter(item) {
                // Check name for valid extension while Galaxy is sniffing the file type
                let ext = item.name.match(/[^.]+$/);
                ext = ext ? ext[0] : '';
                return this.format.includes(item.extension)
                    || (item.extension === '' && this.format.includes(ext))
                    || (item.extension === 'auto' && this.format.includes(ext))
                    || (item.extension === 'data' && this.format.includes(ext));
            },
            onInput(models) {
                this.selected_models = models;
                if (models.length === 0) {
                    this.selection = null;
                } else if (this.type === 'data_collection_input' || models.length > 1) {
                    this.selection = {
                        src: 'new_collection',
                        collection_type: 'list',
                        name: this.$attrs.id, // TODO id is provided by input annotation JSON
                        element_identifiers: models.reduce((acc, model) => {
                            // Convert model to input object
                            let result = model.toInput();

                            // Ensure name is unique
                            let name = model.name;
                            let suff = 1;
                            while (acc.some(i=>i.name === name)) {
                                name = model.name + '_' + suff.toString();
                                ++suff;
                            }
                            result.name = name;

                            acc.push(result);
                            return acc;
                        }, []),
                    };
                } else {
                    this.selection = models[0].toInput();
                }

                let message = this.selection;
                if (!message && this.$attrs.id) {
                    // TODO remove when optional inputs available
                    message = {subinput: this.$attrs.id};
                } else {
                    message.subinput = this.$attrs.id;
                }
                this.$emit('input', message);
            },
            setCustomValidity(message){
                /* Sets a custom validity message for the element. If this message is not the empty string,
                then the element is suffering from a custom validity error, and does not validate.
                 */
                this.validation_message = message;
            },
            checkValidity() {
                for (const model of this.selected_models) {
                    if (model.state !== 'ok') return false;
                }
                if (this.optional && !this.selection) return true;
                return this.validator(this.selection) === '' && this.validation_message === '';
            },
            reportValidity() {
                for (const model of this.selected_models) {
                    if (model.state !== 'ok') {
                        this.validation_message = "Some selected datasets are not ready for analysis. Please wait or make a different selection.";
                        return false;
                    }
                }
                if (this.optional && !this.selection) return true;
                this.validation_message = this.validator(this.selection);
                return this.validation_message === '';
            },
            reset() {
                this.$refs.history_contents.clearSelected();
            },
        },
    }
</script>

<style scoped>

</style>