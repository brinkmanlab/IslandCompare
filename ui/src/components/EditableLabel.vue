<template>
    <span class="editable-label">
        <template v-if="edited_label === null">{{ value }}</template>
        <b-input v-else
                 v-model="edited_label"
                 v-bind:pattern="pattern"
                 v-bind:formatter="value=>value.trim()"
                 v-bind:placeholder="placeholder"
                 autofocus
                 @keypress.enter="stop_edit"
                 @blur="stop_edit"
                 @click.stop
                 ref="editbox"
        ></b-input>
    </span>
</template>

<script>
    export default {
        name: "EditableLabel",
        props: {
            value: {
                type: String,
                required: true,
            },
            placeholder: {
                type: String,
                default: '',
            },
            pattern: {
                type: String,
                default: '.+',
            },
        },
        data: ()=>{return{
            edited_label: null,
        }},
        methods: {
            start_edit() {
                this.edited_label = this.value;
                this.$nextTick(()=>{
                    this.$refs.editbox.focus();
                    this.$refs.editbox.select();
                })
            },
            stop_edit() {
                if (this.edited_label !== null) {
                    this.$emit('update', this.edited_label);
                }
                this.edited_label = null;
            },
        },
    }
</script>

<style scoped>
    .editable-label input {
        font-size: inherit;
        padding: 0;
        height: unset;
    }
</style>