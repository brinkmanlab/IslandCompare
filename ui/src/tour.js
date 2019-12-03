// https://pulsar.gitbooks.io/vue-tour/

export const steps = (tutorial)=>[
    {
        target: '.galaxy-workflow-parameters [order="0"]',
        content: 'Begin by uploading your data to be analysed by dragging and dropping it into the box to the right. Alternatively, click the upload button <i class="icon-file-upload"></i> and select your datasets to upload.',
        offset: -100,
        params: {
            placement: 'left',
        }
    },
    {
        target: '.galaxy-workflow-parameters [order="0"]',
        content: 'Once your datasets have been uploaded, select them by clicking in the box to the right. Hold Ctrl (⌘ for mac) to select multiple. Hold Shift to select a range.',
        offset: -100,
        params: {
            placement: 'left',
        }
    },
    {
        target: '.galaxy-workflow-parameters [order="1"]',
        content: 'Optionally, you can upload a newick formatted tree which will decide the order of alignment in the visualization. If no tree is provided, one will be automatically computed.',
        offset: -100,
        params: {
            placement: 'left',
        }
    },
    {
        target: '.galaxy-workflow-parameters [order="2"]',
        content: 'If you are analysing draft genomes, select a reference genome that is closest to your datasets species. The drafts will be stitched in the order of alignment to the reference. If no reference is selected the draft contigs will be stitched in the order they appear in the uploaded dataset.',
        offset: -100,
        params: {
            placement: 'left',
        }
    },
    {
        target: '.invocation-name',
        content: 'Enter a label to distinguish this analysis from others.',
        offset: -100,
        params: {
            placement: 'bottom-left',
        }
    },
    {
        target: '.invocation-submit',
        content: 'Click submit to run the analysis.',
        offset: -100,
        params: {
            onCreate(data) { //https://popper.js.org/popper-documentation.html#onUpdate
                const target = data.instance.reference;
                target.addEventListener('click', tutorial.nextStep);
            },
            placement: 'bottom-left',
        }
    },
    {
        target: '.analysis-tabs',
        content: 'The pending job will appear in the Recent Jobs tab. Once complete a "Visualize" button will appear in the Job History page along with the option to download the analysis. <br/><a href="/visualize?src=https%3A%2F%2Fislandcompare.pathogenomics.sfu.ca%2Fdemo%2Flisteria_sample_analysis.gff3" target="_self">View Example</a>',
        offset: -100,
        params: {
            placement: 'right',
        }
    },
];

export const callbacks = {
    onStart() {
    },
    onPreviousStep(step) { //eslint-disable-line
    },
    onNextStep(step) { //eslint-disable-line
    },
    onStop() {
    },
};

export default {
    steps, callbacks,
}