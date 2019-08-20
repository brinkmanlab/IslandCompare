<template>
    <!-- HDA, HDCA -->
    <DatasetParameter v-if="['data_input', 'data_collection_input'].includes(type)"
                      v-bind:label="label"
                      v-bind:type="type"
                      v-bind:annotation="annotation"
                      v-bind:historyPromise="historyPromise"
                      v-bind:value="value || []"
                      @input="onInput"
                      @upload="onUpload"
                      ref="param"
    />
    <!-- Reference Genome -->
    <DBKeyParameter v-else-if="type === 'dbkey_input'"
                    v-bind:label="label"
                    v-bind:type="type"
                    v-bind:annotation="annotation"
                    v-bind:value="value"
                    @input="onInput"
                    ref="param"
    />
    <!-- Simple parameter -->
    <SimpleParameter v-else-if="type === 'parameter_input'"
                     v-bind:label="label"
                     v-bind:type="type"
                     v-bind:annotation="annotation"
                     v-bind:value="value"
                     @input="onInput"
                     ref="param"
    />
    <span v-else>Unrecognised input type</span>
</template>

<script>
    import DatasetParameter from "./WorkFlowParameters/DatasetParameter";
    import DBKeyParameter from "./WorkFlowParameters/DBKeyParameter";
    import SimpleParameter from "./WorkFlowParameters/SimpleParameter";

    export default {
        name: "WorkflowParameter",
        components: {
            SimpleParameter,
            DBKeyParameter,
            DatasetParameter,
        },
        props: {
            type: {
                type: String,
                required: true,
            },
            label: {
                type: String,
                required: true,
            },
            uuid: {
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
                //TODO default: History.getRecent(),
            },
            value: {
                type: [String, Number, Array],
                default: '',
            }
        },
        methods: {
            onInput(value) {
                this.$emit('input', { [this.uuid]: value });
            },
            onUpload(upload) {
                // Uploads handled here rather than in DatasetParameter as it may be generalised for tools later.
                for (const file of upload.files) {
                    let ext = file.name.match(/[^.]+$/);
                    if (upload.permitted_file_extensions.length === 0 || (ext && upload.permitted_file_extensions.includes(ext[0]))) {
                        upload.history.fileUpload(file);
                    } else {
                        let tmp_id = file.name + Math.floor(Math.random() * 10 ** 16).toString();
                        upload.history.insert({
                            data: {
                                id: tmp_id,
                                file: file,
                                name: "Incorrect file format: " + file.name,
                                hid: -1,
                                history_id: this.history.id
                            }
                        });
                    }
                }
            },
            setCustomValidity(message){
                /* Sets a custom validity message for the element. If this message is not the empty string,
                then the element is suffering from a custom validity error, and does not validate.
                 */
                this.$refs.param.setCustomValidity(message);
            },
            checkValidity() {
                /* Returns a Boolean that is false if the element is a candidate for constraint validation,
                and it does not satisfy its constraints. In this case, it also fires an invalid event at the element.
                It returns true if the element is not a candidate for constraint validation, or if it satisfies
                its constraints.
                 */
                return this.$refs.param.checkValidity();
            },
            reportValidity() {
                /* Runs the checkValidity() method, and if it returns false (for an invalid input or no pattern
                attribute provided), then it reports to the user that the input is invalid in the same manner
                as if you submitted a form.
                 */
                return this.$refs.param.reportValidity();
            },
            reset() {
                this.$refs.param.reset();
            },
        }
    }
</script>

<style scoped>

</style>