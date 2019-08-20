<template>
    <b-card class="galaxy-workflow-parameter-dataset" v-bind:header="label" v-bind:border-variant="validation_message ? 'danger' : 'default'">
        <HistoryContents ref="history_contents"
                         v-bind:historyPromise="historyPromise"
                         v-bind:value="value"
                         v-bind:permitted_file_extensions="permitted_file_extensions"
                         @input="onInput"
                         @upload="upload=>$emit('upload', upload)"
        />
        <span v-if="mapped">
            <i class="icon-batch-mode"></i>
            This input expects a single selection. Separate jobs will be triggered for each dataset selection.
        </span>
        <b-card-footer v-if="validation_message" footer-text-variant="danger">
            <em>{{validation_message}}</em>
        </b-card-footer>
    </b-card>
</template>

<script>
    import HistoryContents from "../../histories/HistoryContents";
    export default {
        name: "DatasetParameter",
        components: {HistoryContents},
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
            validator: {
                type: Function,
                default(selection) {
                    return selection.length > 0 ? '' : 'Please select at least one dataset';
                }
            }
        },
        data() {return{
            selection: [],
            validation_message: '',
        }},
        computed: {
            mapped() {
                return this.type === 'data_input' && this.selected.length > 1;
            },
            permitted_file_extensions() {
                //Workflow inputs have no type specification. If a json array of strings is present in the input annotation, use those.
                const match = this.annotation.match(/\[[^\]]+]/);
                if (match) {
                    return JSON.parse(match[0]);
                }
                return [];
            },
        },
        methods: {
            onInput(models) {
                this.selection = models.map(model=>model.toInput());
                this.$emit('input', this.selection);
            },
            setCustomValidity(message){
                /* Sets a custom validity message for the element. If this message is not the empty string,
                then the element is suffering from a custom validity error, and does not validate.
                 */
                this.validation_message = message;
            },
            checkValidity() {
                /* Returns a Boolean that is false if the element is a candidate for constraint validation,
                and it does not satisfy its constraints. In this case, it also fires an invalid event at the element.
                It returns true if the element is not a candidate for constraint validation, or if it satisfies
                its constraints.
                 */
                return this.validator(this.selection) === '';
            },
            reportValidity() {
                /* Runs the checkValidity() method, and if it returns false (for an invalid input or no pattern
                attribute provided), then it reports to the user that the input is invalid in the same manner
                as if you submitted a form.
                 */
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