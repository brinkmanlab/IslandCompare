$(document).ready(function(){
    loadGenomesToTable();
    loadJobsToTable();
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
                tablerowbuilder = tablerowbuilder.concat("<td>"+
                    "<input type=\"checkbox\" name=\"jobCheckList\" value="+data[index].id+" />"+"</td>");
                tablerowbuilder = tablerowbuilder.concat("</tr>");
                genomeTable.append(tablerowbuilder);
            }
        }
    })
}

function loadJobsToTable(){
    $.ajax({
        type:"GET",
        url:"getJobs",
        success: function(data){
            jobsTable = $("#jobTable > tbody");
            $("#jobTable tbody > tr").remove(); // Clear all data from UI, this will be reloaded
            for (var index=0;index<data.length;index++){
                tablerowbuilder = "";
                tablerowbuilder = tablerowbuilder.concat("<tr>");
                for (var item in data[index]){
                    tablerowbuilder = tablerowbuilder.concat("<td>"+data[index][item]+"</td>");
                }
                tablerowbuilder = tablerowbuilder.concat("<td>"+"<a href=\'"+"getAlignment/"+
                    data[index].id+"\'"+">Mauve</a></td>");
                tablerowbuilder = tablerowbuilder.concat("</tr>");
                jobsTable.append(tablerowbuilder);
            }
        }
    })
}

$(document).ready(function(){
    //Send Job Request to Server
    $("#genomeListForm").submit(function(){
        $.post("/submitJob",
            $("#genomeListForm").serialize(),
            function(response){
                console.log(response);
            });
        return false;
    });
});