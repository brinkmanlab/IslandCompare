$(document).ready(function(){
    loadGenomesToTable();

    // Listener for Genome Upload form submission
    $("#genomeUploadForm").submit(function(){
        // TODO Write form submission code here, return false on completion
        loadGenomesToTable();
    });
});

// Load Genomes into Status Table in UI
function loadGenomesToTable(){
    $.ajax({
        type:"GET",
        url: "getGenomes",
        success: function(data){
            genomeTable = $("#genomeTable > tbody");
            $("#genomeTable tbody > tr").remove(); // Clear all data from UI, this will be reloaded
            for (var index=0;index<data.length;index++){
                delete data[index]['uploader']; // Remove uploader from output dict
                tablerowbuilder = "";
                tablerowbuilder = tablerowbuilder.concat("<tr>");
                for (var item in data[index]){
                    tablerowbuilder = tablerowbuilder.concat("<td>"+data[index][item]+"</td>");
                }
                tablerowbuilder = tablerowbuilder.concat("</tr>");
                genomeTable.append(tablerowbuilder);
            }
        }
    })
}