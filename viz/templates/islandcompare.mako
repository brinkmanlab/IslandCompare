<html>
<head>
    <script>
        let query = new URLSearchParams(location.search);
        let dataset_id = query.get('dataset_id');
        query.delete('dataset_id');
        query.set('src', '/galaxy/datasets/' + dataset_id + '/display');
        window.location.href = window.location.href.replace(/\/islandcompare\/show.*/, '/islandcompare/static/islandcompare.html?' + query.toString());
    </script>
</head>
</html>