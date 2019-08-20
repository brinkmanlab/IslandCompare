<template>
    <!-- Boolean -->
    <b-form-radio-group v-if="this.type === 'boolean'" buttons
                        class="galaxy-workflow-parameter-simple"
                        v-bind:state="validation_message ? 'invalid' : null"
    >
        <option value="true">Yes</option>
        <option value="false">No</option>
        <b-form-invalid-feedback>{{validation_message}}</b-form-invalid-feedback>
    </b-form-radio-group>
    <!-- Text, Integer, Float, Color -->
    <b-form-group v-else
                  v-bind:label="label"
                  v-bind:description="description"
                  v-bind:invalid-feedback="validation_message"
                  v-bind:state="validation_message ? 'invalid' : null"
                  class="galaxy-workflow-parameter-simple"
    >
        <b-form-input v-bind:type="form_type" v-bind:number="form_type === 'number'" ref="input"/>
    </b-form-group>
</template>

<script>
    export default {
        name: "SimpleParameter",
        props: {
            type: {
                type: String,
                required: true,
            },
            label: {
                type: String,
                required: true,
            },
            description: {
                type: String,
                default: '',
            }
        },
        data() {return{
            validation_message: '',
        }},
        computed: {
            form_type() {
                // Map Galaxy input types to HTML input types
                switch (this.type) {
                    default:
                    case 'boolean':
                    case 'text':
                        return 'text';
                    case 'integer':
                    case 'float':
                        return 'number';
                    case 'color':
                        return 'color';
                }
            },
        },
        methods: {
            setCustomValidity(message){
                /* Sets a custom validity message for the element. If this message is not the empty string,
                then the element is suffering from a custom validity error, and does not validate.
                 */
                this.validation_message = message;
                if (this.type !== 'boolean') this.$refs.input.setCustomValidity(message);
            },
            checkValidity() {
                /* Returns a Boolean that is false if the element is a candidate for constraint validation,
                and it does not satisfy its constraints. In this case, it also fires an invalid event at the element.
                It returns true if the element is not a candidate for constraint validation, or if it satisfies
                its constraints.
                 */
                if (this.type === 'boolean') return true; // If a boolean value has only one valid state, why does it exist?

            },
            reportValidity() {
                /* Runs the checkValidity() method, and if it returns false (for an invalid input or no pattern
                attribute provided), then it reports to the user that the input is invalid in the same manner
                as if you submitted a form.
                 */
                if (this.type !== 'boolean') this.$refs.input.reportValidity();
            },
        },
    }
</script>

<style scoped>

</style>