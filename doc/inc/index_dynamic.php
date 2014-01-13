<?php include('faqdata.php'); ?>
<?php

?>
<style type="text/css">

p.q {
  font-weight:   bold;
}

ul {
  list-style-type: circle;
}

h4 {
  font-size: 14pt;
  color: #5050C0;
  padding-top: 3em;
}

.reveal:hover {
    color:#ccb;
}
.reveal {
    padding-top: 1em;
    cursor: pointer;
    margin: 5px 0;
}

.answer {
    border: 1px dotted grey;
    background-color: #eef;
    padding-left:   1em;
}

</style>

<html>
<head>
<script src="http://code.jquery.com/jquery-1.10.2.js"></script>
<script src="jquery.transit.min.js"></script>

</head>
<body>
<?php
//Output a table of contents with popping up

echo '<div id="faq">';
$sectionCounter = 1;
foreach ($faq as $section) {
  echo '<h4>'.$section['title'].'</h4><br />';
  echo '<ol>';
    foreach($section['content'] as $item) {
    $label = $item['label'];
  ?>
  <li><div id="<?php echo $label ?>" class="reveal" name="#<?php echo $label ?>">
      <?php echo $item['question']; ?>
      </div>
       <div id="<?php echo $label ?>_ans" class="answer">
       <?php echo $item['answer'];?>
      </div>
  </li>
  <?php
  }
  echo "  </ol>";
}




/*
//Output the faq
foreach ($faq as $section) {
  echo '<h4>'.$section['title'].'</h4>';
  foreach($section['content'] as $item) {
    $label=$item['label'];
    echo "<a name=\"$label\" />";
  ?>
  <p class="q"><?php echo $item['question'];?> </p>
  <p class="a"><?php echo $item['answer'];?> </p>


  <?php
  }

}*/

echo '</div>';
?>
<script>
$(".answer").hide();
$(".reveal").click(function() {
    var questionId = "#"+this.id+"_ans";
    $(questionId).toggle("fast");
});

</script>


</body>

</html>