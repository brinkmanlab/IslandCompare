<template>
    <ul>
        <li v-for="genome in genomes"
            v-bind:key="genome.id">
            <span>{{ genome.name }}</span>
        </li>
    </ul>
</template>

<script>
    import * as galaxy from '@/galaxy'
    export default {
        name: "ReferenceGenomes",
        computed: {
            genomes() {
                return galaxy.genomes.Genome.getAll();
            }
        },
        beforeMount() {
            galaxy.genomes.Genome.$fetch()
        }
    }
</script>

<style scoped>

</style>