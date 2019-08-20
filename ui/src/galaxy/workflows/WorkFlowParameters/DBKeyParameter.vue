<template>
    <b-card class="galaxy-workflow-parameter-dbkey" v-bind:border-variant="validation_message ? 'danger' : 'default'">
        <b-card-header header-tag="header">{{label}}</b-card-header>
        <ReferenceGenomes v-bind:value="value" @input="onInput" ref="table"></ReferenceGenomes>
        <b-card-footer v-if="validation_message" footer-text-variant="danger">
            <em>{{validation_message}}</em>
        </b-card-footer>
    </b-card>
</template>

<script>
    import ReferenceGenomes from "../../genomes/ReferenceGenomes";
    export default {
        name: "DBKeyParameter",
        components: {ReferenceGenomes},
        props: {
            label: {
                type: String,
                required: true,
            },
            value: {
                type: String,
                default: '',
            }
        },
        data() {return {
            validation_message: '',
            selection: '',
        }},
        methods: {
            onInput(value) {
                this.selection = value;
                this.$emit('input', value);
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
                return this.selection !== '';
            },
            reportValidity() {
                /* Runs the checkValidity() method, and if it returns false (for an invalid input or no pattern
                attribute provided), then it reports to the user that the input is invalid in the same manner
                as if you submitted a form.
                 */
                if (this.checkValidity()) return true;
                else {
                    this.validation_message = "Please select a reference dataset";
                    return false;
                }
            },
            reset() {
                this.$refs.table.clearSelected();
            },
        }
    }
</script>

<style scoped>

</style>